import gzip
import heapq
import multiprocessing
import queue
import random
import sys  # Needed for the exit function
import time
from heapq import heappop, heappush
from queue import Empty

import pyzstd
from rich import console
from rich.console import group
from rich.panel import Panel
from rich.progress import Progress, track


def read_file(filename, read_start_event, queue, reader_barrier, read_complete_event):
    read_start_event.set()
    sequences = []
    count = 0

    chunk_number = 0
    with gzip.open(filename, 'rt') as f:
        for line in f:
            count += 1
            if count % 4 == 2:
                sequences.append(line.strip())
                
                if len(sequences) == 1000000:
                    reader_barrier.wait()
                    queue.put((chunk_number, sequences.copy()))
                    chunk_number += 1
                    sequences.clear()
                    console.log(f"[bold purple]Reader {multiprocessing.current_process().name}[/bold purple] proceeding to read next chunk", style="yellow")

    if sequences:
        reader_barrier.wait()
        queue.put((chunk_number, sequences.copy()))  # Added the chunk_number metadata here
        console.log(f"[bold purple]Reader {multiprocessing.current_process().name}[/bold purple] completed reading chunk of size: ", len(sequences), style="green")
        sequences.clear()


    # Signal the completion of the entire read for this reader
    read_complete_event.set()
    console.log(f"[bold purple]Reader {multiprocessing.current_process().name}[/bold purple] completed all reads, shutting down.", style="red")



def sorter(input_queues, sorted_chunks_queue, read_complete_events, sort_complete_event, sorter_lock):    
    while True:
        all_queues_empty = all([q.empty() for q in input_queues])
        all_readers_done = all([event.is_set() for event in read_complete_events])

        if all_queues_empty and all_readers_done:
            break

        chunks = []

        with sorter_lock:
            try:
                chunk_tuples = [q.get(timeout = 0.5) for q in input_queues]
                
                # Separate the chunk_number from the data
                chunk_numbers, chunks = zip(*chunk_tuples)
                
                # Check if all chunk_numbers are the same. If not, there's an error.
                if len(set(chunk_numbers)) != 1:
                    console.log(f"[bold white]Sorter {multiprocessing.current_process().name}[/bold white] encountered mismatched chunk_numbers.", style="red")
                    continue

                tuples_to_sort = list(zip(*chunks))
                console.log(f"[bold white]Sorter {multiprocessing.current_process().name}[/bold white] got chunks of sizes: {', '.join(map(str, [len(chunk) for chunk in chunks]))} with chunk_number {chunk_numbers[0]}", style="yellow")
            except queue.Empty:
                continue

        sorted_tuples = sorted(tuples_to_sort)
        sorted_chunks_queue.put(sorted_tuples)
        
        console.log(f"[bold white]Sorter {multiprocessing.current_process().name}[/bold white] completed sorting chunk sizes: {', '.join(map(str, [len(chunk) for chunk in chunks]))}", style="green")

    sort_complete_event.set()
    console.log(f"[bold white]Sorter {multiprocessing.current_process().name}[/bold white] completed all sorting, shutting down.", style="red")

def merger(sorted_chunks_queue, output_queues, sort_complete_events, merger_done_event):
    heap = []  # Local heap to manage chunks based on their sizes
    
    total_sorters = len(sort_complete_events)
    finished_sorters = 0
    
    while finished_sorters < total_sorters or sorted_chunks_queue.qsize() >= 2 or len(heap) >= 2:
        queue_chunks = []
        finished_sorters = sum([event.is_set() for event in sort_complete_events])

        # If there are 2 or more chunks in the queue, merge them first
        if sorted_chunks_queue.qsize() >= 2:
            console.log("A")
            chunks_got = 0
            while not sorted_chunks_queue.empty():
                chunk = sorted_chunks_queue.get()
                queue_chunks.extend(chunk)
                chunks_got+=1

        if queue_chunks:
            console.log("B")
            queue_merged_chunk = sorted(queue_chunks)
            console.log(f"[bold blue]Merger {multiprocessing.current_process().name}[/bold blue] merged {chunks_got} queue chunks from queue.", style="green")
            heappush(heap, (len(queue_merged_chunk), queue_merged_chunk))

        elif len(heap) >= 2:
            console.log("C")
            _, chunk1 = heappop(heap)
            _, chunk2 = heappop(heap)
            merged_chunk = sorted(chunk1 + chunk2)
            console.log(f"[bold blue]Merger {multiprocessing.current_process().name}[/bold blue] merged heap chunks of sizes: {len(chunk1)}, {len(chunk2)}", style="green")
            heappush(heap, (len(merged_chunk), merged_chunk))
            continue

    # Check if any sorters have completed
        # console.log("D")
    if(finished_sorters == total_sorters):
        console.log("D")
        console.log(f"[bold blue]Merger {multiprocessing.current_process().name}[/bold blue] detected all sorters have completed.", style="yellow")
            # When the heap contains only 1 chunk, check the queue one last time
    while not sorted_chunks_queue.empty():
        console.log("E")
        chunk = sorted_chunks_queue.get()
        heappush(heap, (len(chunk), chunk))
        console.log(f"[bold blue]Merger {multiprocessing.current_process().name}[/bold blue] pulled chunk from queue to the heap.", style="green")

    # Merge remaining chunks from the heap
    while len(heap) > 1:
        console.log("F")    
        _, chunk1 = heappop(heap)
        _, chunk2 = heappop(heap)
        merged_chunk = sorted(chunk1 + chunk2)
        console.log(f"[bold blue]Merger {multiprocessing.current_process().name}[/bold blue] merged heap chunks of sizes: {len(chunk1)}, {len(chunk2)}", style="green")
        heappush(heap, (len(merged_chunk), merged_chunk))

    console.log(f"[bold blue]Merger {multiprocessing.current_process().name}[/bold blue] completed merging. Final chunk size: {len(heap[0][1]) if heap else 0}.", style="green")

    # Directly pop from heap, split the data, and send to the output queues
    while heap:
        console.log("G")
        console.log(f"[bold blue]Merger {multiprocessing.current_process().name}[/bold blue] is preparing to send data to writers.", style="yellow")
        _, chunk = heappop(heap)
        separated_lists = list(zip(*chunk))
        console.log(f"[bold blue]Merger {multiprocessing.current_process().name}[/bold blue] separated chunk of size: ", len(chunk), style="green")
        for i, q in enumerate(output_queues):
            q.put(separated_lists[i])
            console.log(f"[bold blue]Merger {multiprocessing.current_process().name}[/bold blue] sent chunk of size: ", len(separated_lists[i]), style="green")
    
    merger_done_event.set()



def write_output(queue, filename, merger_done_event):
    # Wait for the merger to finish
    merger_done_event.wait()

    sequences = []
    while True:
        try:
            chunk = queue.get(timeout = 1)  # Try to get data from the queue with a timeout
            sequences.extend(chunk)
            console.log(f"{filename.replace('.fastq.gz', '.reads.zst')} writer {multiprocessing.current_process().name} got chunk of size: ", len(chunk), style="yellow")
        except Empty:  # If no data is available in the timeout period
            if queue.empty():  # If the queue is empty, break out of the loop
                break

    # Use pyzstd to compress and write to .zst file directly
    with pyzstd.open(filename.replace('.fastq.gz', '.reads.zst'), 'wb') as f:
        for seq in sequences:
            f.write((seq + '\n').encode())
    
    console.log(f"{filename.replace('.fastq.gz', '.reads.zst')} writer {multiprocessing.current_process().name} completed writing, shutting down.", style="red")

if __name__ == "__main__":
    # Configuration
    num_sorters = 5
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