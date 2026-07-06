import sys
import csv
import gzip
import mmh3
import math
from array import array
from collections import defaultdict, deque, Counter
from typing import List, Tuple, Optional, Set, Dict, Generator

# -- PARAMETERS -- #

K = 13      # K-mer length
W = 26      # Window size for winnowing

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

def load_fasta_gz(filepath: str) -> Generator[str, None, None]:
   """Yields reads from .fasta.gz file."""
   try:
      with gzip.open(filepath, 'rt') as f:
         for line in f:
            line = line.strip()
            if line and not line.startswith('>'):
               yield line.upper() 
   except FileNotFoundError:
      print(f"[ERROR] File {filepath} not found.")
   except Exception as e:
      print(f"[ERROR] Something went wrong when reading {filepath}: {e}")

# -- ALGORITHMS -- #

UINT64_MASK = (1 << 64) - 1

def to_uint64(number: int) -> int:
   return number & UINT64_MASK

def hash_kmer(kmer: str) -> int:
   lo, hi = mmh3.hash64(kmer, signed=True)
   return to_uint64(to_uint64(lo) ^ ((to_uint64(hi) << 1) & UINT64_MASK))

def get_minimizers(seq: str) -> Set[int]:
   """Compute minimizer set for a single sequence using the winnowing scheme."""

   # Check if creating k-mer is impossible.
   if len(seq) < K:
      return set()

   # Create list of hashed kmers.
   n_kmers = len(seq) - K + 1
   kmer_hashes = []
   for i in range(n_kmers):
      kmer = seq[i:i+K]
      kmer_hashes.append(hash_kmer(kmer))

   # Check if window is larger than number of k-mers.
   if W > n_kmers:
      return {min(kmer_hashes)}

   minimizers = set()   # Minimal hashes
   dq = deque()         # Indices of hashed kmers

   # Winnowing scheme:
   for i in range(n_kmers):
      # Remove all hashes that are larger than current one.
      while dq and kmer_hashes[i] < kmer_hashes[dq[-1]]:
         dq.pop()

      # Add current index to the last position.
      dq.append(i)

      # Remove all indices outside the window.
      while dq and dq[0] <= i - W:
         dq.popleft()

      # Add new minimizer if the window is about to move or has already moved.
      if i >= W - 1:
         min_hash = dq[0]
         minimizers.add(kmer_hashes[min_hash])

   return minimizers

def calculate_score(sample_set: Set[int], class_array: array) -> float:
   """Calculates the similarity score between a query sample and a known class."""
   if not sample_set or not class_array:
      return 0.0
   inter = sum(1 for m in class_array if m in sample_set)
   return 0.5 * (inter / len(sample_set) + inter / len(class_array))

# ----- MAIN ----- #

def main():
   if len(sys.argv) != 4:
      print("Usage: python3 classifier.py training_data.tsv testing_data.tsv output.tsv")
      sys.exit(1)
   
   training_file, testing_file, output_file = sys.argv[1], sys.argv[2], sys.argv[3]

   print("Loading training metadata...")
   training_datasets = load_tsv_data(training_file, TRAINING_COLUMNS)
   print(f"Loaded {len(training_datasets)} training datasets.")

   # Group files by class to process them together.
   files_by_class = defaultdict(list)
   for file, loc in training_datasets:
      files_by_class[loc].append(file)
   del training_datasets

   class_minimizers: Dict[str, array] = {}

   # Build Class Sketches - Training Phase
   max_n_reads = 0 #? Move to loop
   for loc, files in files_by_class.items():
      print(f"Processing class {loc}...") #!
      loc_counter = Counter()

      # Process all files belonging to this class
      for file in files:
         n_reads = 0
         reads = load_fasta_gz(file)
         print(f"Training file {file}...") #!
         for read in reads:
            if read:
               loc_counter.update(get_minimizers(read))
               n_reads += 1
         max_n_reads = max(max_n_reads, n_reads)

      # Noise filtering, calculate threshold based on read amount
      minimizer_limit = 2 ** math.log10(max(1, max_n_reads // 1000))
      valid_minimizers = [m for m, c in loc_counter.items() if c >= minimizer_limit]
      del loc_counter

      valid_minimizers.sort()
      class_minimizers[loc] = array('Q', valid_minimizers) # Compact array
      del valid_minimizers
      print(f"Class {loc}: kept {len(class_minimizers[loc])} unique minimizers, limit = {minimizer_limit}.") #!
   del files_by_class

   locations = sorted(class_minimizers.keys())
   print(f"Classes found: {locations}.") #!

   print("Loading testing metadata...")
   testing_datasets = load_tsv_data(testing_file, TESTING_COLUMNS)
   print(f"Loaded {len(testing_datasets)} testing datasets.")

   # Classify Test Samples - Testing Phase
   results = []
   for (file, ) in testing_datasets:
      reads = load_fasta_gz(file)
      print(f"Testing file {file}...") #!

      # Extract minimizers from test sample
      sample_minimizers: Set[int] = set()
      for read in reads:
         if read:
            sample_minimizers.update(get_minimizers(read))

      # Calculate similarity scores against all trained classes
      scores = {}
      for loc in locations:
         scores[loc] = calculate_score(sample_minimizers, class_minimizers[loc])
      print(f"Scores for {file}: {scores}.") #!
      results.append((file, scores))

   with open(output_file, 'w', encoding='utf-8') as out:
      header = ["fasta_file"] + locations
      out.write("\t".join(header) + "\n")
      for file, scores in results:
         row = [file] + [f"{scores.get(loc, 0.0):.6f}" for loc in locations]
         out.write("\t".join(row) + "\n")

if __name__ == "__main__":
   main()
