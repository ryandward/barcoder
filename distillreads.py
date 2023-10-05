import multiprocessing
import gzip
import sys

def read_file(filename, queue, read_complete_event):
    sequences = []
    count = 0

    with gzip.open(filename, 'rt') as f:
        for line in f:
            count += 1
            if count % 4 == 2:
                sequences.append(line.strip())
                if len(sequences) == 40000000:
                    queue.put(sequences.copy())
                    sequences.clear()

    if sequences:
        queue.put(sequences)
    read_complete_event.set()

def sorter(input_queues, output_queues, read_complete_events):
    while True:
        # Wait for all readers to have at least one chunk or be done
        while not all([not q.empty() for q in input_queues]) and not all([event.is_set() for event in read_complete_events]):
            pass

        # Check if all readers are done and all queues are empty
        if all([event.is_set() for event in read_complete_events]) and all([q.empty() for q in input_queues]):
            break

        # Get chunks from all queues to create tuples
        chunks = [q.get() for q in input_queues]

        # Create tuples for sorting
        tuples_to_sort = list(zip(*chunks))

        # Sort the tuples
        sorted_tuples = sorted(tuples_to_sort)

        # Split the tuples and put in respective output queues
        for i, queue in enumerate(output_queues):
            queue.put([tup[i] for tup in sorted_tuples])
    # Signal completion to writers
    for q in output_queues:
        q.put(None)


def write_output(queue, filename, read_complete_event):
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
    filenames = sys.argv[1:]
    manager = multiprocessing.Manager()

    queues = [manager.Queue() for _ in filenames]
    read_complete_events = [multiprocessing.Event() for _ in filenames]
    output_queues = [manager.Queue() for _ in filenames]

    readers = [multiprocessing.Process(target=read_file, args=(filename, queue, event)) for filename, queue, event in zip(filenames, queues, read_complete_events)]
    sorter_process = multiprocessing.Process(target=sorter, args=(queues, output_queues, read_complete_events))
    writers = [multiprocessing.Process(target=write_output, args=(queue, filename, event)) for queue, filename, event in zip(output_queues, filenames, read_complete_events)]

    for r in readers:
        r.start()
    
    sorter_process.start()

    for w in writers:
        w.start()

    for r in readers:
        r.join()

    sorter_process.join()

    for w in writers:
        w.join()
