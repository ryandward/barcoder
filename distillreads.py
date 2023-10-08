import gzip
import heapq
import multiprocessing as mp
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

#import defaultdict
from collections import defaultdict

def read_file(filename, chunk_complete_sem, queue, reader_barrier, read_complete_event):
    sequences = []
    count = 0
    chunk_number = 0
    reader_name = f"[bold purple]Reader {mp.current_process().name}[/bold purple]"
    console.log(f"{reader_name} beginning to read: {filename}", style="yellow")

    with gzip.open(filename, 'rt') as f:
        for line in f:
            count += 1
            if count % 4 == 2:
                sequences.append(line.strip())
                
                if len(sequences) == 2**19:
                    queue.put((chunk_number, sequences.copy()))
                    console.log(f"{reader_name} read {((count//4)):,} sequences, chunk: {chunk_number:,}", style="white")
                    chunk_number += 1
                    sequences.clear()


    if sequences:
        queue.put((chunk_number, sequences.copy()))
        console.log(f"{reader_name} finished at {((count//4)):,} sequenced, chunk: {chunk_number:,}", style="bold green")
        sequences.clear()


    # Signal the completion of the entire read for this reader
    read_complete_event.set()
    console.log(f"{reader_name} finished.", style="red")



def dispatch(dispatch_queues, chunk_complete_sems, sorter_queue, read_complete_events, dispatch_done_event):    
    dispatch_name = f"[bold magenta]Dispatch {mp.current_process().name}[/bold magenta]"
    chunks_dict = defaultdict(list)

    while True:
        all_readers_done = all([event.is_set() for event in read_complete_events])

        # Exit condition
        if all([q.empty() for q in dispatch_queues]) and all_readers_done:
            break

        # Fetch available chunks from each reader
        for idx, q in enumerate(dispatch_queues):
    
            chunk_number, chunk_data = q.get()
            # Use idx to label the origin of the chunk
            chunks_dict[chunk_number].append((idx, chunk_data))

        # Check chunks_dict for any complete sets of chunks
        complete_chunk_numbers = [cn for cn, data in chunks_dict.items() if len(data) == len(dispatch_queues)]
        
        for cn in complete_chunk_numbers:
            # Strip out the index label and sort based on the index to ensure correct order
            sorted_data = sorted(chunks_dict[cn], key=lambda x: x[0])
            stripped_data = [item[1] for item in sorted_data]
            tuples_to_sort = list(zip(*stripped_data))
            sorter_queue.put(tuples_to_sort)

            sequence_lengths = [len(seq) for seq in stripped_data]
            console.log(f"{dispatch_name} dispatched a chunk to sorters: {cn} with sequence lengths: {'; '.join(['{:,}'.format(length) for length in sequence_lengths])}", style="green")
            del chunks_dict[cn]
            
        time.sleep(0.2)

    for event in read_complete_events: event.wait()
    dispatch_done_event.set()
    console.log(f"{dispatch_name} finished.", style="red")




def sorter(sorter_queue, chunk_complete_sems, sorted_chunks_queue, read_complete_events, sort_complete_event, sorter_lock):    
    sorter_name = f"[bold white]Sorter {mp.current_process().name}[/bold white]"

    while True:
        tuples_to_sort = []
        sorter_queue_empty = sorter_queue.empty()
        dispatch_done = dispatch_done_event.is_set()

        if sorter_queue_empty and dispatch_done:
            break

        with sorter_lock:

            try:
                tuples_to_sort = sorter_queue.get(timeout = 0.1)

            except queue.Empty:
                pass
                # console.log(f"{sorter_name} did not retrieve chunks, trying again.", style="yellow")
                
            finally:
                pass
                # Release semaphores in the finally block to ensure they are always released

        # Moved outside the lock block
        if tuples_to_sort:  # Checking to ensure there's data to process
            console.log(f"{sorter_name} starting a chunk: {len(tuples_to_sort[0]):,} x {len(tuples_to_sort):,}", style="yellow")
            sorted_tuples = sorted(tuples_to_sort)
            sorted_chunks_queue.put(sorted_tuples)
            console.log(f"{sorter_name} finished a chunk: {len(tuples_to_sort[0]):,} x {len(tuples_to_sort):,}", style="green")

    for event in read_complete_events: event.wait()
    dispatch_done_event.wait()
    sort_complete_event.set()

    console.log(f"{sorter_name} finished.", style="red")


def merger(sorted_chunks_queue, output_queue_filenames, sort_complete_events, merger_done_event, writer_sems):
    heap = []
    total_sorters = len(sort_complete_events)
    finished_sorters = 0
    merger_name = f"[bold blue]Merger {mp.current_process().name}[/bold blue]"

    while finished_sorters < total_sorters or sorted_chunks_queue.qsize() >= 2 or len(heap) >= 2:
        queue_chunks = []
        finished_sorters = sum([event.is_set() for event in sort_complete_events])

        if sorted_chunks_queue.qsize() >= 2:
            while not sorted_chunks_queue.empty():
                chunk = sorted_chunks_queue.get()
                queue_chunks.extend(chunk)

        if queue_chunks:
            queue_merged_chunk = sorted(queue_chunks)
            console.log(f"{merger_name} merged {len(queue_chunks):,} sequences from queue.", style="green")
            heappush(heap, (len(queue_merged_chunk), queue_merged_chunk))
        elif len(heap) >= 2:
            _, chunk1 = heappop(heap)
            _, chunk2 = heappop(heap)
            merged_chunk = sorted(chunk1 + chunk2)
            console.log(f"{merger_name} merged heap sequences of sizes: {len(chunk1):,} and {len(chunk2):,}.", style="bold green")
            heappush(heap, (len(merged_chunk), merged_chunk))

    if finished_sorters == total_sorters:
        console.log(f"{merger_name} confirmed all sorters have completed.", style="yellow")

    for event in read_complete_events: event.wait()
    dispatch_done_event.wait()
    for event in sort_complete_events: event.wait()

    while not sorted_chunks_queue.empty():
        chunk = sorted_chunks_queue.get()
        heappush(heap, (len(chunk), chunk))
        console.log(f"{merger_name} moved {len(chunk):,} sequences from queue to heap.", style="bold yellow")

    while len(heap) > 1:
        _, chunk1 = heappop(heap)
        _, chunk2 = heappop(heap)
        merged_chunk = sorted(chunk1 + chunk2)
        console.log(f"{merger_name} merged remaining heap sequences of sizes: {len(chunk1):,} and {len(chunk2):,}.", style="bold green")
        heappush(heap, (len(merged_chunk), merged_chunk))

    console.log(f"{merger_name} final merged chunk size: {len(heap[0][1]):,}.", style="bold green")

    merger_done_event.set()
    console.log(f"{merger_name} contacting writers...", style="bold white")    
    # Dispatch to writers
    while heap:
        _, chunk = heappop(heap)
        separated_lists = list(zip(*chunk))
        for i, (q, filename) in enumerate(output_queue_filenames):
            q.put(separated_lists[i])
            console.log(f"{merger_name} dispatched {len(separated_lists[i]):,} sequences to writer-{i+1} for file {filename}.", style="green")
            writer_sems[i].release()

    console.log(f"{merger_name} finished.", style="red")

def write_output(queue, filename, merger_done_event, writer_sem):
    merger_done_event.wait()

    sequences = []
    output_filename = filename.replace('.fastq.gz', '.reads.zst')
    writer_name = f"[bold cyan]Writer {mp.current_process().name}[/bold cyan]"

    while True:
        try:
            writer_sem.acquire()
            chunk = queue.get(timeout=1)
            sequences.extend(chunk)
            console.log(f"{writer_name} retrieved {len(chunk):,} sequences.", style="yellow")

            with pyzstd.open(output_filename, 'wb') as f:
                for seq in sequences:
                    f.write((seq + '\n').encode())
            console.log(f"{writer_name} finished writing {len(sequences):,} to ({output_filename}). Shutting down.", style="red")
        except Empty:
            if queue.empty():
                break
        finally:
            try:
                writer_sem.release()
            except ValueError:
                pass  # Semaphore was already released, no action needed

def write_output(queue, filename, merger_done_event, writer_sem):
    merger_done_event.wait()

    sequences = []
    output_filename = filename.replace('.fastq.gz', '.reads.zst')
    writer_name = f"[bold cyan]Writer {mp.current_process().name}[/bold cyan]"

    while True:
        try:
            writer_sem.acquire()
            chunk = queue.get(timeout=1)
            sequences.extend(chunk)
            console.log(f"{writer_name} retrieved a chunk of {len(chunk):,} sequences.", style="yellow")

            with pyzstd.open(output_filename, 'wb') as f:
                for seq in sequences:
                    f.write((seq + '\n').encode())
            console.log(f"{writer_name} finished writing {len(sequences):,} sequences to ({output_filename}). Shutting down.", style="red")
        except Empty:
            if queue.empty():
                break
        finally:
            try:
                writer_sem.release()
            except ValueError:
                console.log(f"{writer_name} already released semaphore, synchronizing with merger.", style="yellow")




if __name__ == "__main__":
    # Configuration
    num_sorters = 12
    filenames = sys.argv[1:]

    from rich.console import Console
    from rich.progress import Progress, SpinnerColumn, TimeElapsedColumn

    console = Console(record=True, stderr=True)

    # Set up mp resources
    with mp.Manager() as manager:

        chunk_complete_sems = [manager.Semaphore(0) for _ in filenames]
        
        dispatch_queues = [manager.Queue() for _ in filenames]
        sorter_queue = manager.Queue()

        read_complete_events = [manager.Event() for _ in filenames]
        dispatch_done_event = manager.Event()

        sort_start_events = [manager.Event() for _ in range(num_sorters)]
        sorter_lock = manager.Lock()
        sorted_chunks_queue = manager.Queue()
        sort_complete_events = [manager.Event() for _ in range(num_sorters)]

        merger_lock = manager.Lock()
        merger_done_event = manager.Event()
        
        output_queue_filenames = [(manager.Queue(), filename) for filename in filenames]

        writer_sems = [manager.Semaphore(0) for _ in range(len(output_queue_filenames))]


        # Create processes
        reader_barrier = mp.Barrier(len(filenames))
        readers = [mp.Process(target=read_file, args=(filename, chunk_complete_sem, queue, reader_barrier, event)) for filename, chunk_complete_sem, queue, event in zip(filenames, chunk_complete_sems, dispatch_queues, read_complete_events)]
        dispatch_process = mp.Process(target=dispatch, args=(dispatch_queues, chunk_complete_sems, sorter_queue, read_complete_events, dispatch_done_event))
        sorters = [mp.Process(target=sorter, args=(sorter_queue, chunk_complete_sems, sorted_chunks_queue, read_complete_events, sort_complete_event, sorter_lock)) for sort_complete_event in sort_complete_events]
        merger_process = mp.Process(target=merger, args=(sorted_chunks_queue, output_queue_filenames, sort_complete_events, merger_done_event, writer_sems))
        writers = [mp.Process(target=write_output, args=(queue, filename, merger_done_event, writer_sem)) for (queue, filename), writer_sem in zip(output_queue_filenames, writer_sems)]



        # Start processes
        for r in readers:
            r.start()
        dispatch_process.start()
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