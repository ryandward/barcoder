import gzip
import heapq
import multiprocessing
import queue
import random
import sys
import time
from queue import Empty
from rich.progress import track
import pyzstd
from rich import console
from rich.console import group
from rich.panel import Panel
from rich.progress import Progress

def read_file(filename, read_start_event, queue, reader_barrier, read_complete_event):
    read_start_event.set()
    sequences = []
    count = 0

    with gzip.open(filename, 'rt') as f:
        for line in f:
            count += 1
            if count % 4 == 2:
                sequences.append(line.strip())
                if len(sequences) == 1048576*2:
                    queue.put(sequences.copy())
                    console.log(f"Reader {multiprocessing.current_process().name} completed reading chunk of size: ", len(sequences), style="bold green")
                    sequences.clear()
                    reader_barrier.wait()
                    console.log(f"Reader {multiprocessing.current_process().name} proceeding to read next chunk", style="bold yellow")

    if sequences:
        queue.put(sequences)

    # Signal the completion of the entire read for this reader
    read_complete_event.set()
    console.log(f"Reader {multiprocessing.current_process().name} completed all reads, shutting down.", style="bold red")

def sorter(input_queues, sorted_chunks_queue, read_complete_events, sort_complete_event, sorter_lock):
    # console.log(f"Sorter {multiprocessing.current_process().name} started.")
    
    while True:
        all_queues_empty = all([q.empty() for q in input_queues])
        all_readers_done = all([event.is_set() for event in read_complete_events])

        # If all queues are empty and all readers are done, break out of the loop
        if all_queues_empty and all_readers_done:
            break

        chunks = []
        with sorter_lock:
            try:
                # Try to get chunks from all queues with a timeout
                chunks = [q.get(timeout = 1) for q in input_queues]
                tuples_to_sort = list(zip(*chunks))
                console.log(f"Sorter {multiprocessing.current_process().name} got chunks of sizes: {', '.join(map(str, [len(chunk) for chunk in chunks]))}", style="bold yellow")
                # Create tuples for sorting
            except queue.Empty:
                continue
   
        # Sort the tuples
        sorted_tuples = sorted(tuples_to_sort)
        # Send the sorted tuples to the merger
        sorted_chunks_queue.put(sorted_tuples)
        
        console.log(f"Sorter {multiprocessing.current_process().name} completed sorting chunk sizes: {', '.join(map(str, [len(chunk) for chunk in chunks]))}", style="bold green")

    sort_complete_event.set()
    console.log(f"Sorter {multiprocessing.current_process().name} completed all sorting, shutting down.", style="bold red")

def k_way_merge(chunks):
    merged = []
    heap = []

    for idx, chunk in enumerate(chunks):
        if chunk:  # If chunk is not empty
            heapq.heappush(heap, (chunk[0], idx, 0))
    
    while heap:
        val, chunk_idx, element_idx = heapq.heappop(heap)
        merged.append(val)
        
        if element_idx + 1 < len(chunks[chunk_idx]):
            next_val = chunks[chunk_idx][element_idx + 1]
            heapq.heappush(heap, (next_val, chunk_idx, element_idx + 1))
    
    return merged

def merger(sorted_chunks_queue, output_queues, sort_complete_events, merger_done_event):
    total_sorters = len(sort_complete_events)
    finished_sorters = 0
    merged_data = []  # This will hold the continuously merged data

    while finished_sorters < total_sorters:
        try:
            chunk = sorted_chunks_queue.get(timeout = 1)
            if chunk:
                console.log(f"Merger {multiprocessing.current_process().name} got chunk of size: {len(chunk)}", style="bold yellow")
                
                # Merge the current chunk with the existing merged_data
                merged_data = sorted(merged_data + chunk)
                console.log(f"Merger {multiprocessing.current_process().name} merged data size: {len(merged_data)}", style="bold green")
                
        except Empty:
            finished_sorters = sum([event.is_set() for event in sort_complete_events])

    # Splitting tuples into separate lists for final output
    separated_lists = list(zip(*merged_data))
    for i, q in enumerate(output_queues):
        console.log(f"Merger {multiprocessing.current_process().name} sending chunk of size: {len(separated_lists[i])} to writer queue.", style="bold yellow")
        q.put(separated_lists[i])
        console.log(f"Merger {multiprocessing.current_process().name} sent chunk of size: {len(separated_lists[i])} to writer queue.", style="bold green")
    merger_done_event.set()
    console.log(f"Merger {multiprocessing.current_process().name} completed all merging, shutting down.", style="bold red")


def write_output(queue, filename, merger_done_event):
    # Wait for the merger to finish
    merger_done_event.wait()

    sequences = []
    while True:
        try:
            chunk = queue.get(timeout = 1)  # Try to get data from the queue with a timeout
            sequences.extend(chunk)
            console.log(f"{filename.replace('.fastq.gz', '.reads.zst')} writer {multiprocessing.current_process().name} got chunk of size: ", len(chunk), style="bold yellow")
        except Empty:  # If no data is available in the timeout period
            if queue.empty():  # If the queue is empty, break out of the loop
                break

    # Use pyzstd to compress and write to .zst file directly
    with pyzstd.open(filename.replace('.fastq.gz', '.reads.zst'), 'wb') as f:
        for seq in sequences:
            f.write((seq + '\n').encode())
    
    console.log(f"{filename.replace('.fastq.gz', '.reads.zst')} writer {multiprocessing.current_process().name} completed writing, shutting down.", style="bold red")

if __name__ == "__main__":
    # Configuration
    num_sorters = 12
    filenames = sys.argv[1:]

    from rich.console import Console
    from rich.progress import Progress, SpinnerColumn, TimeElapsedColumn

    console = Console(record=True, stderr=True)

    # Set up multiprocessing resources
    with multiprocessing.Manager() as manager:

        read_start_events = [manager.Event() for _ in filenames]
        reader_queues = [manager.Queue() for _ in filenames]
        read_complete_events = [manager.Event() for _ in filenames]

        sort_start_events = [manager.Event() for _ in range(num_sorters)]
        sorter_lock = manager.Lock()
        sorted_chunks_queue = manager.Queue()
        sort_complete_events = [manager.Event() for _ in range(num_sorters)]

        merger_lock = manager.Lock()
        merger_done_event = manager.Event()
        output_queues = [manager.Queue() for _ in filenames]

        # Create processes
        reader_barrier = multiprocessing.Barrier(len(filenames))
        readers = [multiprocessing.Process(target=read_file, args=(filename, read_start_event, queue, reader_barrier, event)) for filename, read_start_event, queue, event in zip(filenames, read_start_events, reader_queues, read_complete_events)]
        sorters = [multiprocessing.Process(target=sorter, args=(reader_queues, sorted_chunks_queue, read_complete_events, sort_complete_events[i], sorter_lock)) for i in range(num_sorters)]
        merger_process = multiprocessing.Process(target=merger, args=(sorted_chunks_queue, output_queues, sort_complete_events, merger_done_event))
        writers = [multiprocessing.Process(target=write_output, args=(queue, filename, merger_done_event)) for queue, filename in zip(output_queues, filenames)]



        # Start processes
        for r in readers:
            r.start()
        for s in sorters:
            s.start()
        merger_process.start()
        for w in writers:
            w.start()
            
        # # Wait for processes to complete
        # for r in readers:
        #     r.join()
        # for s in sorters:
        #     s.join()
        # merger_process.join()
        for w in writers:
            w.join()