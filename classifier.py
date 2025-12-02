import time
import sys
import csv
import gzip
from typing import List, Tuple, Dict
from Bio import SeqIO

# ----- TIME ----- #

def now() -> float:
   return time.time()

def since(t0: float) -> str:
   return "%.4f" % (time.time() - t0)

# ----- FILE ----- #

# Important columns in .tsv files
TRAINING_COLUMNS = [0, 1]
TESTING_COLUMNS = [0]

def load_tsv_data(filepath: str, required_indices: List[int]) -> List[Tuple[str, str]]:
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
               row_data = (row[index] for index in required_indices)
               data.append(row_data)

   except FileNotFoundError:
      print(f"[ERROR] File {filepath} not found.")
      sys.exit(1)

   return data

def load_fasta_gz(filepath: str):
   pass
   # TODO Load thingy

# ----- CODE ----- #

# TODO Other thingys

# ----- MAIN ----- #

def main():
   if len(sys.argv) != 4:
      print("Usage: python3 classifier.py training_data.tsv testing_data.tsv output.tsv")
      sys.exit(1)
      
   training_file, testing_file, output_file = sys.argv[1], sys.argv[2], sys.argv[3]

   print("Loading training metadata...")
   training_datasets = load_tsv_data(training_file, TRAINING_COLUMNS)
   print(f"Loaded {len(training_datasets)} training datasets.")
   if training_datasets:
      print(f"Example dataset: {training_datasets[0]}")

   print("Loading testing metadata...")
   testing_datasets = load_tsv_data(testing_file, TESTING_COLUMNS)
   print(f"Loaded {len(testing_datasets)} testing datasets.")
   if testing_datasets:
      print(f"Example dataset: {testing_datasets[0]}")

   # TODO Main thingy

   print(f"Saving output: {output_file}")
   with open(output_file, 'w', encoding='utf-8') as outfile:
      outfile.write("fasta_file\tpredicted_class\n")
      for fasta_file in testing_datasets:
         outfile.write(f"{fasta_file}\tPLACEHOLDER_CLASS\n")

if __name__ == "__main__":
   main()