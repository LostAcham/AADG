import time
import sys
import csv
import gzip
from typing import List, Tuple, Optional
import mmh3
import numpy as np
from collections import defaultdict

# ----- TIME ----- #

def now() -> float:
   return time.time()

def since(t0: float) -> str:
   return "%.4f" % (time.time() - t0)

# ----- FILE ----- #

# Important columns in .tsv files
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

# ----- CODE ----- #

STEP = 1                         # Kmer step
K = 4                            # Kmer length
SEEDS = [i for i in range(200)]  # Number of hashes

def generate_hashed_kmers(seq: str, seed: int) -> List[str]:
   """Generates hashed kmers from a sequence."""
   return set(mmh3.hash(seq[pos:pos + K], seed) for pos in range(0, len(seq) - K + 1, STEP))

def generate_min_hashes(seq: str) -> List[int]:
   """Generates set of min hashes."""
   return [min(generate_hashed_kmers(seq, seed)) for seed in SEEDS]

def calculate_centroid(signatures: List[List[int]]) -> Optional[List[int]]:
   if not signatures:
      return None
   return np.min(signatures, axis=0)

def minhash_similarity(sig1: List[int], sig2: List[int]) -> float:
   """Calculates similarities between two signatures (based on minHashes)."""
   if sig1 is None or sig2 is None:
      return 0.0
   equal = 0
   for hash1, hash2 in zip(sig1, sig2):
      if hash1 == hash2:
         equal += 1
   return equal / len(sig1)

# ----- MAIN ----- #

def main():
   if len(sys.argv) != 4:
      print("Usage: python3 classifier.py training_data.tsv testing_data.tsv output.tsv")
      sys.exit(1)
      
   training_file, testing_file, output_file = sys.argv[1], sys.argv[2], sys.argv[3]

   print("Loading training metadata...")
   training_datasets = load_tsv_data(training_file, TRAINING_COLUMNS)
   print(f"Loaded {len(training_datasets)} training datasets.")

   # TODO Training
   sketches = defaultdict(list)
   for (file, loc) in training_datasets:
      loaded_reads = load_fasta_gz(file)
      print(f"Loaded {len(loaded_reads)} training reads.")
      for loaded_read in loaded_reads:
         if not loaded_read:
            continue
         sketches[loc].append(generate_min_hashes(loaded_read))
   
   locations = sorted(sketches.keys())
   centroids = defaultdict(list)
   for loc in locations:
      centroids[loc] = calculate_centroid(sketches[loc])

   print("Loading testing metadata...")
   testing_datasets = load_tsv_data(testing_file, TESTING_COLUMNS)
   print(f"Loaded {len(testing_datasets)} testing datasets.")

   # TODO Testing
   results = []
   for (file, ) in testing_datasets:
      loaded_reads = load_fasta_gz(file)
      print(f"Loaded {len(loaded_reads)} testing reads.")
      if not loaded_reads:
         results.append((file, {loc: 0.0 for loc in locations}))
         continue
      signature = []
      for loaded_read in loaded_reads:
         signature.append(generate_min_hashes(loaded_read))
      signature = list(map(min, zip(*signature)))
      sims = [minhash_similarity(signature, centroids[loc]) for loc in locations]
      scores = {c: p for c, p in zip(locations, sims)}
      print(scores)
      results.append((file, scores))

   # TODO Actual writing after testing is done
   print(f"Saving output: {output_file}")
   with open(output_file, 'w', encoding='utf-8') as outfile:
         header = ["fasta_file"] + locations
         outfile.write("\t".join(header) + "\n")
         for (f, scores) in results:
            row = [f] + [f"{scores.get(c,0.0):.6f}" for c in locations]
            outfile.write("\t".join(row) + "\n")

if __name__ == "__main__":
   main()