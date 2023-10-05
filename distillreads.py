import multiprocessing
import gzip
import sys
import random
import time

def read_file(filename, queue, barrier, read_complete_event):
    sequences = []
    count = 0

    with gzip.open(filename, 'rt') as f:
        for line in f:
            count += 1
            if count % 4 == 2:
                sequences.append(line.strip())
                if len(sequences) == 65536*4:
                    queue.put(sequences.copy())
                    print(f"Reader {multiprocessing.current_process().name} completed reading chunk of size: ", len(sequences))
                    sequences.clear()
                    barrier.wait()  # Wait for all other readers here
                    print(f"Reader {multiprocessing.current_process().name} proceeding to read next chunk")

    if sequences:
        queue.put(sequences)

    # Signal the completion of the entire read for this reader
    read_complete_event.set()
    print(f"Reader {multiprocessing.current_process().name} completed all reads, shutting down.")


def sorter(input_queues, sorted_chunks_queue, read_complete_events, sort_complete_event, sorter_lock):
    print(f"Sorter {multiprocessing.current_process().name} started.")
    while True:
        all_queues_empty = all([q.empty() for q in input_queues])
        all_readers_done = all([event.is_set() for event in read_complete_events])

        if all_readers_done:
            print(f"Sorter {multiprocessing.current_process().name} detected all readers are done.")
        if all_queues_empty:
            print(f"Sorter {multiprocessing.current_process().name} detected all queues are empty.")

        # Check if all readers are done and all queues are empty
        if all([event.is_set() for event in read_complete_events]) and all([q.empty() for q in input_queues]):
            break
        
        chunks = []
        # Only proceed if all the queues are not empty
        with sorter_lock:  # Acquire the lock to ensure exclusive access
            if all([not q.empty() for q in input_queues]):
                # Get chunks from all queues to create tuples
                chunks = [q.get() for q in input_queues]
                print(f'Sorter {multiprocessing.current_process().name} got chunks of size: ', [len(chunk) for chunk in chunks]) 
                 
            # Create tuples for sorting
            tuples_to_sort = list(zip(*chunks))

            # Sort the tuples
            sorted_tuples = sorted(tuples_to_sort)
            print('Sorted chunks of size: ', len(sorted_tuples))

            # Send the sorted tuples to the merger
            sorted_chunks_queue.put(sorted_tuples)
            print(f"Sorter {multiprocessing.current_process().name} completed sorting chunk size: ", len(sorted_tuples))

    print(f"Sorter {multiprocessing.current_process().name} reader events: {[event.is_set() for event in read_complete_events]}")
    print(f"Sorter {multiprocessing.current_process().name} queues: {[q.empty() for q in input_queues]}")
    print(f"Sorter {multiprocessing.current_process().name} completed all sorting, shutting down.")
    sort_complete_event.set()
    sorted_chunks_queue.put(None)  # Signal completion to merger
   
    print(f"Sorter {multiprocessing.current_process().name} completed all sorting and set the complete event.")






import heapq

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
        # Wait for all sorters to complete their tasks
    while not all([event.is_set() for event in sort_complete_events]):
        pass

    all_chunks = []

    # Counter for the number of end signals (None values) received
    end_signals_received = 0

    while end_signals_received < len(sort_complete_events):
        chunk = sorted_chunks_queue.get()
        print(f"end_signals_received: {end_signals_received}", len(sort_complete_events))
        
        # Check for end signal
        if chunk is None:
            end_signals_received += 1
        elif len(chunk) > 0:
            print('Merger got chunk of size: ', len(chunk))
            all_chunks.append(chunk)
            print(f'Merger has {len(all_chunks)} chunks in total')

    print('Merger got all chunks, starting k-way merge of size: ', len(all_chunks))
    merged_list = k_way_merge(all_chunks)
    print('K-way merged a list of size: ', len(merged_list))

    # Splitting tuples into separate lists
    separated_lists = list(zip(*merged_list))
    for i, queue in enumerate(output_queues):
        queue.put(separated_lists[i])

    # Signal completion to writers
    for q in output_queues:
        q.put(None)
    
    merger_done_event.set()



def write_output(queue, filename, merger_done_event):
    # Wait for the merger to finish
    merger_done_event.wait()
    
    sequences = []
    while True:
        chunk = queue.get()
        if chunk is None:  # Check for end signal
            break
        sequences.extend(chunk)

    with open(filename.replace('.fastq.gz', '.reads'), 'w') as f:
        for seq in sequences:
            f.write(seq + '\n')


if __name__ == "__main__":
    num_sorters = 100  # Or any other number based on available resources

    filenames = sys.argv[1:]
    manager = multiprocessing.Manager()
    sorter_lock = manager.Lock()
    merger_done_event = manager.Event()


    queues = [manager.Queue() for _ in filenames]
    chunk_complete_events = [multiprocessing.Event() for _ in filenames]
    read_complete_events = [multiprocessing.Event() for _ in filenames]
    sort_complete_events = [multiprocessing.Event() for _ in range(num_sorters)]



    output_queues = [manager.Queue() for _ in filenames]
    sorted_chunks_queue = manager.Queue()


    merger_process = multiprocessing.Process(target=merger, args=(sorted_chunks_queue, output_queues, sort_complete_events, merger_done_event))
    barrier = multiprocessing.Barrier(len(filenames))  # One barrier point for each reader
    readers = [multiprocessing.Process(target=read_file, args=(filename, queue, barrier, event)) for filename, queue, event in zip(filenames, queues, read_complete_events)]
    sorters = [multiprocessing.Process(target=sorter, args=(queues, sorted_chunks_queue, read_complete_events, sort_complete_events[i], sorter_lock)) for i in range(num_sorters)]
    writers = [multiprocessing.Process(target=write_output, args=(queue, filename, merger_done_event)) for queue, filename in zip(output_queues, filenames)]

    for r in readers:
        r.start()

    merger_process.start()
    for s in sorters:
        s.start()

    for w in writers:
        w.start()

    for r in readers:
        r.join()

    for s in sorters:
        s.join()

    merger_process.join()

    for w in writers:
        w.join()
