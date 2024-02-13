import gzip
import heapq
import multiprocessing as mp
import os
import queue
import random
import sys
import time
from collections import defaultdict
from heapq import heappop, heappush
from multiprocessing import cpu_count
from queue import Empty

import psutil
import pyzstd
from rich.console import Console
from rich.console import group
from rich.panel import Panel
from rich.progress import Progress, track

console = Console()


def zstd_compress(sequences):
    data = "\n".join(sequences)
    compressed_data = pyzstd.compress(data.encode("utf-8"))
    return compressed_data


def zstd_decompress(compressed_data):
    decompressed_data = pyzstd.decompress(compressed_data).decode("utf-8")
    return decompressed_data.split("\n")


def check_memory_usage(percentage=True):
    memory_info = psutil.virtual_memory()
    if percentage:
        return memory_info.percent
    else:
        return memory_info.used / (1024**3)  # Return usage in GB


def check_cpu_usage():
    return psutil.cpu_percent(interval=1)


def read_file(filename, dispatch_queues, read_done_event):
    sequences = []
    count = 0
    chunk_number = 0
    reader_name = f"[bold purple]Reader {mp.current_process().name}[/bold purple]"
    console.log(f"{reader_name} beginning to read: {filename}", style="yellow")

    # Check file extension to determine how to open
    open_func = gzip.open if filename.endswith(".gz") else open

    with open_func(filename, "rt") as f:
        for line in f:
            count += 1
            if count % 4 == 2:
                sequences.append(line.strip())
                if len(sequences) == 2**20:
                    compressed_chunk = zstd_compress(sequences.copy())
                    dispatch_queues.put((chunk_number, compressed_chunk))
                    console.log(
                        f"{reader_name} read {((count//4)):,} sequences, chunk: {chunk_number:,}",
                        style="white",
                    )
                    chunk_number += 1
                    sequences.clear()

    if sequences:
        compressed_chunk = zstd_compress(sequences.copy())
        dispatch_queues.put((chunk_number, compressed_chunk))
        console.log(
            f"{reader_name} read {((count//4)):,} sequences, chunk: {chunk_number:,}",
            style="white",
        )
        chunk_number += 1
        sequences.clear()

    # Signal the completion of the entire read for this reader
    read_done_event.set()
    if read_done_event.is_set():
        console.log(f"{reader_name} finished.", style="green")
    else:
        console.log(f"{reader_name} failed to finish.", style="bold red")


def dispatch(dispatch_queues, sorter_queue, read_done_events, dispatch_done_event):
    dispatch_name = f"[bold magenta]Dispatch {mp.current_process().name}[/bold magenta]"
    chunks_dict = defaultdict(list)
    finished_queues = set()  # a set to store indexes of finished queues

    while True:
        all_readers_done = all([event.is_set() for event in read_done_events])

        # Exit condition
        if len(finished_queues) == len(dispatch_queues) and sorter_queue.empty():
            console.log(
                f"{dispatch_name} confirmed all readers and sorter have finished.",
                style="yellow",
            )
            break

        for idx, q in enumerate(dispatch_queues):
            # if this queue's reader is done and the queue is empty, add it to the finished set
            if read_done_events[idx].is_set() and q.empty():
                finished_queues.add(idx)

            # if this queue is not marked as finished, try to get data from it
            if idx not in finished_queues:
                try:
                    chunk_number, compressed_chunk_data = (
                        q.get_nowait()
                    )  # non-blocking get
                    chunks_dict[chunk_number].append((idx, compressed_chunk_data))
                except queue.Empty:
                    pass  # this means the queue is currently empty, but the reader might add more to it later

        # Check chunks_dict for any done sets of chunks
        done_chunk_numbers = [
            chunk_number
            for chunk_number, data in chunks_dict.items()
            if len(data) == len(dispatch_queues)
        ]
        # print(done_chunk_numbers)

        for chunk_number in done_chunk_numbers:
            # Sort based on the index to ensure correct order and then strip out the index label
            ordered_chunks = sorted(chunks_dict[chunk_number], key=lambda x: x[0])
            _, chunk_data = list(zip(*ordered_chunks))
            sorter_queue.put(
                chunk_data
            )  # Sending the entire ordered data, sorters will handle decompression

            console.log(
                f"{dispatch_name} dispatched a chunk to sorters: {chunk_number}.",
                style="green",
            )
            del chunks_dict[chunk_number]

    for event in read_done_events:
        event.wait()
    dispatch_done_event.set()
    if dispatch_done_event.is_set():
        console.log(f"{dispatch_name} finished.", style="green")
    else:
        console.log(f"{dispatch_name} failed to finish.", style="bold red")


def sorter(sorter_queue, merger_queue, read_done_events, sort_done_event, sorter_lock):
    sorter_name = f"[bold white]Sorter {mp.current_process().name}[/bold white]"

    while True:
        dispatch_data = []
        sorter_queue_empty = sorter_queue.empty()
        dispatch_done = dispatch_done_event.is_set()

        if sorter_queue_empty and dispatch_done:
            break

        with sorter_lock:
            try:
                dispatch_data = sorter_queue.get(timeout=0.1)
            except queue.Empty:
                pass

        if dispatch_data:
            sorted_data = sorted(
                list(zip(*[zstd_decompress(entry) for entry in dispatch_data]))
            )
            chunk_size = len(sorted_data)

            # peek at the first element of the first 5 sequence to get the size of the elements to make sure they are sorted
            element_sizes = [len(sequence) for sequence in sorted_data[0]]
            merger_queue.put(sorted_data)
            console.log(
                f"{sorter_name} sorted a chunk: {chunk_size:,} x {element_sizes}",
                style="green",
            )

    for event in read_done_events:
        event.wait()
    dispatch_done_event.wait()
    sort_done_event.set()

    console.log(f"{sorter_name} finished.", style="red")


def merger(
    merger_queue,
    send_ends,
    sort_done_events,
    merger_done_event,
    output_filenames,
    writer_sems,
):

    heap = []
    total_sorters = len(sort_done_events)
    finished_sorters = 0
    merger_name = f"[bold blue]Merger {mp.current_process().name}[/bold blue]"

    def merge_from_heap(heap, merger_name):
        _, chunk1 = heappop(heap)
        _, chunk2 = heappop(heap)
        merged_chunk = sorted(chunk1 + chunk2)
        console.log(
            f"{merger_name} merged heap sequences of sizes: {len(chunk1):,} and {len(chunk2):,}.",
            style="bold green",
        )
        heappush(heap, (len(merged_chunk), merged_chunk))

    while (
        finished_sorters < total_sorters or merger_queue.qsize() >= 2 or len(heap) >= 2
    ):
        queue_chunks = []

        finished_sorters = sum([event.is_set() for event in sort_done_events])

        # Handle the chunks in the merger queue
        if merger_queue.qsize() >= 2:
            while not merger_queue.empty():
                chunk = merger_queue.get()
                queue_chunks.extend(chunk)

        # Merge the sequences from the queue
        if queue_chunks:
            queue_merged_chunk = sorted(queue_chunks)
            console.log(
                f"{merger_name} merged {len(queue_chunks):,} sequences from queue.",
                style="green",
            )
            heappush(heap, (len(queue_merged_chunk), queue_merged_chunk))

        # If there are still enough chunks in the heap, merge again
        if len(heap) >= 2:
            merge_from_heap(heap, merger_name)

    for event in read_done_events:
        event.wait()
    dispatch_done_event.wait()
    for event in sort_done_events:
        event.wait()

    while not merger_queue.empty():
        chunk = merger_queue.get()
        heappush(heap, (len(chunk), chunk))
        console.log(
            f"{merger_name} moved {len(chunk):,} sequences from queue to heap.",
            style="bold yellow",
        )

    while len(heap) > 1:
        _, chunk1 = heappop(heap)
        _, chunk2 = heappop(heap)
        merged_chunk = sorted(chunk1 + chunk2)
        console.log(
            f"{merger_name} merged remaining heap sequences of sizes: {len(chunk1):,} and {len(chunk2):,}.",
            style="bold green",
        )
        heappush(heap, (len(merged_chunk), merged_chunk))

    console.log(
        f"{merger_name} final merged chunk size: {len(heap[0][1]):,}.",
        style="bold green",
    )

    console.log(f"{merger_name} contacting writers...", style="bold white")
    # Dispatch to writers

    while heap:
        _, chunk = heappop(heap)
        separated_lists = list(zip(*chunk))

        for i, (send_end, filename) in enumerate(zip(send_ends, output_filenames)):
            # writer_sems[i].release()
            send_end.send(separated_lists[i])  # Send through the pipe
            console.log(
                f"{merger_name} dispatched {len(separated_lists[i]):,} sequences to writer-{i+1} for file {filename}.",
                style="green",
            )

    merger_done_event.set()

    for send_end in send_ends:
        send_end.close()

    # for i, (send_end, filename) in enumerate(zip(send_ends, output_filenames)):
    #         writer_sems[i].acquire()

    console.log(f"{merger_name} finished.", style="red")


def writer_process(recv_end, output_filename, writer_sem):
    writer_name = f"[bold cyan]Writer {mp.current_process().name}[/bold cyan]"

    with pyzstd.open(output_filename, "wb") as f:
        while True:
            # If the merger is done and there's no data in the pipe, break out of the loop
            if merger_done_event.is_set() and not recv_end.poll():
                break

            # Wait for data to be available
            if recv_end.poll():
                # Acquire the semaphore if there's data
                # writer_sem.acquire()

                # Receive the data
                sequences = recv_end.recv()

                # Check for termination condition
                if sequences is None:
                    recv_end.close()
                    break

                # Write the data
                for seq in sequences:
                    f.write((seq + "\n").encode())

                console.log(
                    f"{writer_name} wrote {len(sequences):,} sequences to {output_filename}.",
                    style="green",
                )

    console.log(f"{writer_name} finished.", style="red")


if __name__ == "__main__":
    # Configuration
    num_sorters = cpu_count() // 2
    filenames = sys.argv[1:]

    def get_output_filename(filename):
        if filename.endswith(".fastq.gz"):
            return filename.replace(".fastq.gz", ".reads.zst")
        elif filename.endswith(".fastq"):
            return filename.replace(".fastq", ".reads.zst")
        return filename + ".reads.zst"  # default case

    output_filenames = [get_output_filename(filename) for filename in filenames]

    from rich.console import Console
    from rich.progress import Progress, SpinnerColumn, TimeElapsedColumn

    console = Console(record=True, stderr=True)

    # Set up mp resources
    with mp.Manager() as manager:

        dispatch_queues = [manager.Queue() for _ in filenames]
        sorter_queue = manager.Queue()
        read_done_events = [manager.Event() for _ in filenames]
        dispatch_done_event = manager.Event()
        sort_start_events = [manager.Event() for _ in range(num_sorters)]
        sorter_lock = manager.Lock()
        merger_queue = manager.Queue()
        sort_done_events = [manager.Event() for _ in range(num_sorters)]
        merger_lock = manager.Lock()
        merger_done_event = manager.Event()
        output_queue_filenames = [(manager.Queue(), filename) for filename in filenames]
        writer_sems = [manager.Semaphore(0) for _ in range(len(output_queue_filenames))]

        pipes = [mp.Pipe() for _ in filenames]
        merger_send_ends = [send for send, _ in pipes]
        writer_recv_ends = [recv for _, recv in pipes]

        # Create processes
        readers = [
            mp.Process(target=read_file, args=(filename, queue, event))
            for filename, queue, event in zip(
                filenames, dispatch_queues, read_done_events
            )
        ]
        dispatch_process = mp.Process(
            target=dispatch,
            args=(dispatch_queues, sorter_queue, read_done_events, dispatch_done_event),
        )
        sorters = [
            mp.Process(
                target=sorter,
                args=(
                    sorter_queue,
                    merger_queue,
                    read_done_events,
                    sort_done_event,
                    sorter_lock,
                ),
            )
            for sort_done_event in sort_done_events
        ]
        merger_process = mp.Process(
            target=merger,
            args=(
                merger_queue,
                merger_send_ends,
                sort_done_events,
                merger_done_event,
                output_filenames,
                writer_sems,
            ),
        )
        writers = [
            mp.Process(
                target=writer_process, args=(recv_end, output_filename, writer_sem)
            )
            for recv_end, output_filename, writer_sem in zip(
                writer_recv_ends, output_filenames, writer_sems
            )
        ]

        # Start processes
        for r in readers:
            r.start()
        dispatch_process.start()
        for s in sorters:
            s.start()
        merger_process.start()
        for w in writers:
            w.start()

        # Wait for processes to done
        for r in readers:
            r.join()
        dispatch_process.join()
        for s in sorters:
            s.join()
        merger_process.join()
        for w in writers:
            w.join()
