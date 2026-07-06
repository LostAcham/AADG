#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Made by Adam Kaźmierczak & Marcel Kołakowski
"""

import sys
from typing import List, Tuple
from array import array
from Bio import SeqIO

# ----------------------------------------------- PARAMETERS ----------------------------------------------- #

K = 12 # size of seeds
STEP = 3 # interval between seeds
BIN_SIZE = 250 # size of groups in clustering
CLUSTER_AMOUNT = 3 # quantity of best clusters choosen for edit distance

# -------------------------------------------------- FILE -------------------------------------------------- #

def read_fasta_reference(filename: str) -> str:
	"""Loads reference data from FASTA file."""
	records = list(SeqIO.parse(filename, "fasta"))
	if not records:
		raise ValueError(f"File {filename} is empty.")
	if len(records) > 1:
		print(f"[ALERT] File {filename} has more than once sequence – using only first sequence.")
	return str(records[0].seq)

def read_fasta_reads(filename: str) -> List[Tuple[str, str]]:
	"""Loads list of reads (id, sequence) from FASTA file."""
	reads = []
	for record in SeqIO.parse(filename, "fasta"):
		reads.append((record.id, str(record.seq)))
	if not reads:
		raise ValueError(f"File {filename} is empty.")
	return reads

def write_output(filename: str, mappings: List[Tuple[str, int, int]]) -> None:
	"""Writes output to file in format: read_id \t start \t end \n """
	with open(filename, "w") as f:
		for read_id, start, end in mappings:
			f.write(f"{read_id}\t{start}\t{end}\n")

# ------------------------------------------- KARKKAINEN-SANDERS ------------------------------------------- #

def radixpass(a, b, r, n, k) :
   c = array("I", [0] * (k + 1))
   for i in range(n) :
      c[r[a[i]]] += 1

   somme = 0
   for i in range(k + 1):
      freq, c[i] = c[i], somme
      somme += freq

   for i in range(n) :
      b[c[r[a[i]]]] = a[i]
      c[r[a[i]]] += 1

def suffix_array(s):
   alphabet = [None] + sorted(set(s))
   k = len(alphabet)
   n = len(s)
   t = dict((c, i) for i,c in enumerate(alphabet))
   SA = array('I', [0] * (n + 3))
   kark_sand(array('I', [t[c] for c in s] + [0] * 3), SA, n, k)
   return SA[:n]

def kark_sand(s, SA, n, Km):
   n0  = (n+2) // 3
   n1  = (n+1) // 3
   n2  = n // 3
   n02 = n0 + n2

   SA12 = array('I', [0] * (n02 + 3))
   SA0  = array('I', [0] * n0)

   s12 = [i for i in range(n + (n0 - n1)) if i % 3] 
   s12.extend([0] * 3)
   s12 = array('I', s12)

   radixpass(s12, SA12, s[2:], n02, Km)
   radixpass(SA12, s12, s[1:], n02, Km)
   radixpass(s12, SA12, s, n02, Km)

   name = 0
   c0, c1, c2 = -1, -1, -1
   for i in range(n02) :
      if s[SA12[i]] != c0 or s[SA12[i] + 1] != c1 or s[SA12[i] + 2] != c2 :
         name += 1
         c0 = s[SA12[i]]
         c1 = s[SA12[i] + 1]
         c2 = s[SA12[i] + 2]
      if SA12[i] % 3 == 1 :
         s12[SA12[i] // 3] = name
      else :
         s12[SA12[i] // 3 + n0] = name

   if name < n02 :
      kark_sand(s12, SA12, n02, name + 1)
      for i in range(n02) :
         s12[SA12[i]] = i + 1
   else :
      for i in range(n02) :
         SA12[s12[i] - 1] = i

   s0 = array('I',[SA12[i] * 3 for i in range(n02) if SA12[i] < n0])
   radixpass(s0, SA0, s, n0, Km)

   p = j = k = 0
   t = n0 - n1
   while k < n :
      i = SA12[t] * 3 + 1 if SA12[t] < n0 else (SA12[t] - n0) * 3 + 2
      j = SA0[p] if p < n0 else 0

      if SA12[t] < n0 :
         test = (s12[SA12[t] + n0] <= s12[j // 3]) if(s[i] == s[j]) else (s[i] < s[j])
      elif(s[i] == s[j]) :
         test = s12[SA12[t] - n0 + 1] <= s12[j // 3 + n0] if(s[i + 1] == s[j + 1]) else s[i + 1] < s[j + 1]
      else :
         test = s[i] < s[j]

      if(test) :
         SA[k] = i
         t += 1
         if t == n02 :
            k += 1
            while p < n0 :
               SA[k] = SA0[p]
               p += 1
               k += 1

      else : 
         SA[k] = j
         p += 1
         if p == n0 :
            k += 1
            while t < n02 :
               SA[k] = (SA12[t] * 3) + 1 if SA12[t] < n0 else ((SA12[t] - n0) * 3) + 2
               t += 1
               k += 1
      k += 1

# --------------------------------------------- BWT & FM INDEX --------------------------------------------- #

def build_bwt(text: str, SA: array) -> Tuple[str, int]:
   """Build Burrows-Wheeler Transform from text via suffix array."""
   bw = []
   dollarRow = None
   for i in SA:
      if i == 0:
         dollarRow = len(bw)
         bw.append('$')
      else:
         bw.append(text[i - 1])
   return ''.join(bw), dollarRow

class FmCheckpoints(object):
   """Manages rank checkpoints and handles rank queries, which are O(1) time, 
   with the checkpoints taking O(m) space, where m is length of text."""
   
   def __init__(self, bw, cpIval=8):
      """Scan BWT, creating periodic checkpoints as we go."""
      self.cps = {}        # checkpoints
      self.cpIval = cpIval # spacing between checkpoints
      tally = {}           # tally so far
      # Create an entry in tally dictionary and checkpoint map for
      # each distinct character in text
      for c in bw:
         if c not in tally:
            tally[c] = 0
            self.cps[c] = []
      # Now build the checkpoints
      for i, c in enumerate(bw):
         tally[c] += 1 # up to *and including*
         if i % cpIval == 0:
            for c in tally.keys():
               self.cps[c].append(tally[c])

   def rank(self, bw, c, row):
      """Return # c's there are in bw up to and including row."""
      if row < 0 or c not in self.cps:
         return 0
      i, nocc = row, 0
      # Always walk to left (up) when calculating rank
      while (i % self.cpIval) != 0:
         if bw[i] == c:
            nocc += 1
         i -= 1
      return self.cps[c][i // self.cpIval] + nocc

class FmIndex():
   """O(m) size FM Index, where checkpoints and suffix array samples are spaced O(1) elements apart. 
   Queries like count() and range() are O(n) where n is the length of the query. Finding all k occurrences 
   of a length-n query string takes O(n + k) time."""

   @staticmethod
   def downsampleSuffixArray(sa, n=4):
      """Take only the suffix-array entries for every nth suffix. Keep suffixes at offsets 0, n, 2n... 
      with respect to the text. Return map from the rows to their suffix-array values."""
      ssa = {}
      for i, suf in enumerate(sa):
         if suf % n == 0:
            ssa[i] = suf
      return ssa

   def __init__(self, t, cpIval=8, ssaIval=8):
      if not t.endswith('$'):
         t += '$'
      sa = suffix_array(t)
      self.bwt, self.dollarRow = build_bwt(t, sa) # Get BWT string and offset of $ within it
      self.ssa = self.downsampleSuffixArray(sa, ssaIval) # Get downsampled suffix array
      self.slen = len(self.bwt)
      self.cps = FmCheckpoints(self.bwt, cpIval) # Make rank checkpoints

      # Calculate # occurrences of each character
      tots = {}
      for c in self.bwt:
         tots[c] = tots.get(c, 0) + 1

      # Calculate concise representation of first column
      self.first = {}
      totc = 0
      for c, count in sorted(tots.items()):
         self.first[c] = totc
         totc += count

   def count(self, c):
      """Return number of occurrences of characters < c."""
      if c not in self.first: # c does not occur in text
         for cc in sorted(self.first.keys()):
            if c < cc:
               return self.first[cc]
         return list(self.first.values())[-1]
      else:
         return self.first[c]

   def range(self, seed: str):
      """Return range of BWM rows having seed as a prefix."""
      l, r = 0, self.slen - 1 # closed (inclusive) interval
      for i in range(len(seed) - 1, -1, -1): # from right to left
         c = seed[i]
         l = self.cps.rank(self.bwt, c, l - 1) + self.count(c)
         r = self.cps.rank(self.bwt, c, r) + self.count(c) - 1
         if r < l:
            break
      return l, r + 1

   def resolve(self, row):
      """Given BWM row, return its offset w/r/t T."""
      def stepLeft(row):
         """Step left according to character in given BWT row."""
         c = self.bwt[row]
         return self.cps.rank(self.bwt, c, row - 1) + self.count(c)
      steps = 0
      while row not in self.ssa:
         row = stepLeft(row)
         steps += 1
      return self.ssa[row] + steps

   def hasSubstring(self, seed: str):
      """Return true if and only if seed is substring of indexed text."""
      l, r = self.range(seed)
      return r > l

   def hasSuffix(self, seed: str):
      """Return true if and only if seed is suffix of indexed text."""
      l, r = self.range(seed)
      off = self.resolve(l)
      return r > l and off + len(seed) == self.slen - 1

   def occurrences(self, seed: str):
      """Return offsets for all occurrences of seed, in no particular order."""
      l, r = self.range(seed)
      return [self.resolve(i) for i in range(l, r)]

# -------------------------------------------------- SEED -------------------------------------------------- #

def generate_seeds(read: str) -> List[Tuple[int, str]]:
   """Generates seeds from a read."""
   seeds = []
   for pos in range(0, len(read) - K + 1, STEP):
      seeds.append((pos, read[pos:pos + K]))
   return seeds

def fm_exact(fm: FmIndex, seed: str) -> List[int]:
   """Finds exact match of seed in FM Index."""
   l, r = fm.range(seed)
   if r <= l: # no matches
      return []
   return fm.occurrences(seed)

def cluster_positions(predicted_starts: List[int]) -> List[List[int]]:
   """Sorts predicted starts of seeds into clustered groups."""
   buckets = {}
   for (pred_start) in predicted_starts:
      b = pred_start // BIN_SIZE
      buckets.setdefault(b, []).append(pred_start)

   clusters = []
   for b, vals in buckets.items():
      clusters.append(vals)

   # sort clusters by support descending
   clusters.sort(key=lambda x: len(x), reverse=True)
   return clusters

def find_candidate_windows(read: str, fm: FmIndex) -> List[List[int]]:
   """Generates best clusters of seeds for later computations."""
   seeds = generate_seeds(read)
   predicted_starts = []

   for pos, seed in seeds:
      hits = fm_exact(fm, seed)
      for refpos in hits:
         predicted_starts.append(refpos - pos)

   clusters = cluster_positions(predicted_starts)
   return clusters[:CLUSTER_AMOUNT]

# ------------------------------------------ BANDED EDIT DISTANCE ------------------------------------------ #

def banded_edit_distance(read: str, ref_window: str, ref_window_start: int, band: int, max_cost: int) -> Tuple[int, int, int] | None:
   """Align full read to any substring of reference window using banded edit distance."""
   m = len(read)
   n = len(ref_window)

   prev_jlow = 0 # Lowest index of previous row still in-bound
   prev_jhigh = min(n, band) # Highest index of previous row still in-bound
   prev_row = [0] * (prev_jhigh - prev_jlow + 1) # Previous row (true_index = prev_jlow + index)
   dirs: List[List[int]] = [] # Directions for traceback: 0 - diagonal, 1 - up, 2 - left

   # For every letter in read
   for i in range(1, m + 1):
      j_low = max(0, i - band)
      j_high = min(n, i + band)
      width = j_high - j_low + 1
      curr_row = [max_cost + 1] * width
      curr_dirs = [0] * width

      # For every viable letter in text window
      for j in range(j_low, j_high + 1):
         row_index = j - j_low
         best_cost = max_cost + 1
         best_dir = 0

         # Diagonal (match/mismatch)
         if j - 1 >= prev_jlow:
            best_cost = prev_row[j - 1 - prev_jlow] + (0 if read[i - 1] == ref_window[j - 1] else 1)

         # Up (deletion)
         if j <= prev_jhigh:
            curr_val = prev_row[j - prev_jlow] + 1
            if curr_val < best_cost:
               best_cost = curr_val
               best_dir = 1

         # Left (insertion)
         if j - 1 >= j_low:
            curr_val = curr_row[row_index - 1] + 1
            if curr_val < best_cost:
               best_cost = curr_val
               best_dir = 2

         # Saving the best value
         curr_row[row_index] = best_cost
         curr_dirs[row_index] = best_dir

      # The lowest cost is higher that max cost available
      if min(curr_row) > max_cost:
         return None

      # Saving the row
      dirs.append(curr_dirs)
      prev_row = curr_row
      prev_jlow = j_low
      prev_jhigh = j_high

   # Setting up for traceback
   last_row = prev_row
   j_rel = min(range(len(last_row)), key=lambda x: last_row[x])
   j_end = prev_jlow + j_rel
   cost = last_row[j_rel]
   i = m
   j = j_end

   # Tracebacking until start
   while i > 0:
      row_index = j - max(0, i - band)
      direction = dirs[i - 1][row_index]

      # Diagonal
      if direction == 0:
         i -= 1
         j -= 1
      # Up
      elif direction == 1:
         i -= 1
      # Left
      else:
         j -= 1

   # Converting to absolute coords:
   ref_start = ref_window_start + max(0, j)
   ref_end = ref_window_start + j_end
   return int(cost), int(ref_start), int(ref_end)

# ---------------------------------------------- MAPPING PART ---------------------------------------------- #

def map_read_to_reference(read: str, reference: str, fm: FmIndex) -> Tuple[int, int] | None:
   """Map a single read to the reference data."""
   # Parameters
   READ_LEN = len(read)
   MAX_COST = READ_LEN // 5 - 1 # max edit distance allowed
   BAND_SIZE = MAX_COST + 2 * BIN_SIZE # half-band width in dp table
   WINDOW_SIZE = READ_LEN + BAND_SIZE # extension around predicted region

   # Generating clustering windows
   clusters = find_candidate_windows(read, fm)
   if not clusters:
      return None
   best = None

   # Running banded edit distance algorithm for all clusters
   for cluster in clusters:
      start = max(0, (min(cluster) - min(cluster) % BIN_SIZE) - BIN_SIZE)
      end = min(len(reference), start + WINDOW_SIZE)
      ref_window = reference[start:end]

      result = banded_edit_distance(read, ref_window, start, BAND_SIZE, MAX_COST)

      if result is None:
         continue

      cost, rstart, rend = result

      # Saving the best result
      if (best is None) or (cost < best[0]):
         best = (cost, rstart, rend)

   if best is None:
      return None

   _, rstart, rend = best
   return (rstart, rend)

def map_reads(reference: str, reads: List[Tuple[str, str]], fm: FmIndex) -> List[Tuple[str, int, int]]:
   """Maps every read to reference data."""
   results = []
   for read_id, read in reads:
      mapping = map_read_to_reference(read, reference, fm)
      if mapping:
         start, end = mapping
         results.append((read_id, start, end))
   return results

# -------------------------------------------------- MAIN -------------------------------------------------- #

def main() -> None:
   if len(sys.argv) != 4:
      print("Usage: python3 mapper.py reference_file.fasta reads_file.fasta output_file.txt")
      sys.exit(1)

   ref_file, reads_file, output_file = sys.argv[1], sys.argv[2], sys.argv[3]

   print("[INFO] Loading data...")
   reference = read_fasta_reference(ref_file)
   reads = read_fasta_reads(reads_file)
   print(f"[INFO] Reference length: {len(reference):,} bp.")
   print(f"[INFO] Read quantity: {len(reads)}.")

   print("[INFO] Pre-processing...")
   fm = FmIndex(reference)

   print("[INFO] Mapping...")
   mappings = map_reads(reference, reads, fm)
   print(f"[INFO] Mapped {len(mappings)} / {len(reads)} reads.")

   print("[INFO] Writing results...")
   write_output(output_file, mappings)
   print(f"[INFO] Results written to {output_file}.")

if __name__ == "__main__":
	main()