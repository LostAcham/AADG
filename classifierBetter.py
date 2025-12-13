import sys
import csv
import gzip
import mmh3
import argparse
from typing import List, Tuple, Optional, Set, Dict
from collections import defaultdict, deque
import time

# -- PARAMETERS -- #
K = 13      # K-mer length
W = 11      # Window size for winnowing

# ----- TIME ----- #

def now() -> float:
   return time.time()

def since(t0: float) -> str:
   return "%.4f" % (time.time() - t0)

# ----- FILE ----- #
TRAINING_COLUMNS = [0, 1]
TESTING_COLUMNS = [0]

def load_tsv_data(filepath: str, required_indices: List[int]) -> List[Tuple[str, Optional[str]]]:
   """Loads required data from .tsv file. Ignores header and not needed data."""
   data = []
   try:
      with open(filepath, 'r', encoding='utf-8') as f:
         reader = csv.reader(f, delimiter='\t')

         try: # Skip header
            next(reader)
         except StopIteration:
            print(f"[Error] File {filepath} is empty.")
            sys.exit(1)

         for row in reader:
            if row:
               row_data = tuple(row[index] for index in required_indices)
               data.append(row_data)

   except FileNotFoundError:
      print(f"[ERROR] File {filepath} not found.")
      sys.exit(1)

   return data

def load_fasta_gz(filepath: str) -> List[str]:
   """Loads reads from .fasta.gz file."""
   reads = []
   try:
      with gzip.open(filepath, 'rt') as f:
         for line in f:
            line = line.strip()
            if line and not line.startswith('>'):
               reads.append(line.upper()) 
   except FileNotFoundError:
      print(f"[ERROR] File {filepath} not found.")
      return []
   except Exception as e:
      print(f"[ERROR] Something went wrong when reading {filepath}: {e}")
      return []
      
   return reads

# -- ALGORITHMS -- #
UINT64_MASK = (1 << 64) - 1

def to_uint64(number: int) -> int:
   return number & UINT64_MASK

def hash_kmer(kmer: str) -> int:
   lo, hi = mmh3.hash64(kmer, signed=True)
   return to_uint64(to_uint64(lo) ^ ((to_uint64(hi) << 1) & UINT64_MASK))

def minimizers_for_sequence(seq: str) -> Set[int]:
   """
   Compute minimizer set for a single sequence using the winnowing scheme.
   Returns a set of 64-bit integer hashes (minimizers).
   """

   #! Checking if creating k-mer isn't possible.
   if len(seq) < K:
      return set()

   # Create list of hashed canonicalized kmers.
   nr_of_kmers = len(seq) - K + 1
   kmer_hashes = []
   for i in range(nr_of_kmers):
      kmer = seq[i:i+K]
      kmer_hashes.append(hash_kmer(kmer))

   # Window is larger than number of k-mers.
   if W > nr_of_kmers:
      return set(min(kmer_hashes))

   # Use set to store min_hashes and deque to store indexes of hashed kmers.
   minimizers = set()
   dq = deque()

   for i in range(nr_of_kmers):
      # Remove all kmer hashes that are larger than current one.
      while dq and kmer_hashes[i] < kmer_hashes[dq[-1]]:
         dq.pop()

      # Add current index to the last position.
      dq.append(i)

      # Remove all indexes outside the window.
      while dq and dq[0] <= i - W:
         dq.popleft()

      # Add new minimizer if the window is about to move or has already moved.
      if i >= W - 1:
         min_hash = dq[0]
         minimizers.add(kmer_hashes[min_hash])

   return minimizers

def jaccard_similarity(set_a: Set[int], set_b: Set[int]) -> float:
   if not set_a or not set_b:
      return 0.0
   inter = len(set_a & set_b)
   union = len(set_a | set_b)
   return inter / union

# ----- MAIN ----- #
def main():
   if len(sys.argv) != 4:
      print("Usage: python3 classifier.py training_data.tsv testing_data.tsv output.tsv")
      sys.exit(1)
   
   training_file, testing_file, output_file = sys.argv[1], sys.argv[2], sys.argv[3]

   print("Loading training metadata...")
   training_datasets = load_tsv_data(training_file, TRAINING_COLUMNS)
   print(f"Loaded {len(training_datasets)} training datasets.")

   class_minimizers: Dict[str, Set[int]] = defaultdict(set)

   for file, loc in training_datasets:
      reads = load_fasta_gz(file)
      print(f"Training file {file}: {len(reads)} reads.")
      sample_minimizers: Set[int] = set()
      for read in reads:
         if read:
            sample_minimizers.update(minimizers_for_sequence(read))
      class_minimizers[loc].update(sample_minimizers)

   locations = sorted(class_minimizers.keys())
   print(f"Classes found: {locations}.")

   print("Loading testing metadata...")
   testing_datasets = load_tsv_data(testing_file, TESTING_COLUMNS)
   print(f"Loaded {len(testing_datasets)} testing datasets.")

   results = []
   for (file,) in testing_datasets:
      reads = load_fasta_gz(file)
      print(f"Testing file {file}: {len(reads)} reads.")
      sample_minimizers: Set[int] = set()
      for read in reads:
         if not read:
               continue
         sample_minimizers.update(minimizers_for_sequence(read))

      scores = {}
      for loc in locations:
         scores[loc] = jaccard_similarity(sample_minimizers, class_minimizers[loc])
      print(f"Scores for {file}: {scores}.")
      results.append((file, scores))

   # Write results to output file.
   with open(output_file, 'w', encoding='utf-8') as out:
      header = ["fasta_file"] + locations
      out.write("\t".join(header) + "\n")
      for file, scores in results:
         row = [file] + [f"{scores.get(loc, 0.0):.6f}" for loc in locations]
         out.write("\t".join(row) + "\n")

if __name__ == "__main__":
   main()
