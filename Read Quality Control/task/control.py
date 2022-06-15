import numpy as np
import gzip


class Readfilter:
    def __init__(self, file_path):
        self.lengths = None
        self.n_repeated_reads = None
        self.file_path = file_path
        self.ncontents_filtered = []
        self.ncontents = []
        self.gccontents = []
        self.lines = []
        self.reads = []
        self.get_input(file_path)
        self.extract_reads()
        self.find_repeated_reads()
        self.get_values()

    def get_input(self, file_path):
        with gzip.open(file_path, "r") as file:
            self.lines = [str(line).strip("b'").rstrip('\\n') for line in file]

    def find_repeated_reads(self):
        unique, counts = np.unique(np.array(self.reads), return_counts=True)
        self.n_repeated_reads = sum([count - 1 for count in counts])

    def extract_reads(self):
        self.reads = self.lines[1::4]

    def get_values(self):
        self.lengths = [len(read) for read in self.reads]
        self.gccontents = [round((read.count("G") + read.count("C")) / len(read) * 100, 2) for read in self.reads]
        self.ncontents = [(read.count("N") / len(read)) * 100 for read in self.reads]
        self.ncontents_filtered = list(filter(lambda x: x != 0, self.ncontents))

    def print_results(self):
        print(f"\nReads in the file = {len(self.reads):}")
        print(f"Reads sequence average length = {round(np.mean(self.lengths))}")
        print(f"\nRepeats = {self.n_repeated_reads}")
        print(f"Reads with Ns = {len(self.ncontents_filtered)}")
        print(f"\nGC content average = {round(sum(self.gccontents) / len(self.gccontents), 2)}%")
        print(f"Ns per read sequence = {round(sum(self.ncontents) / len(self.ncontents), 2)}%")


if __name__ == "__main__":
    filepaths = [input(), input(), input()]
    results = []
    for filepath in filepaths:
        results.append(Readfilter(filepath))
    repeats = [readfilter.n_repeated_reads for readfilter in results]
    ncontent = [round(sum(readfilter.ncontents) / len(readfilter.ncontents), 2) for readfilter in results]
    nreads = [len(readfilter.reads) for readfilter in results]
    best = nreads.index(min(nreads))
    results[best].print_results()
