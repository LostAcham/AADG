#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# example mapping program

from collections import defaultdict
import sys
import numpy
from Bio import SeqIO
from sys import argv
import time

start_time = time.time()

# On the presentation tell why only these numbers are possible
POSSIBLE_KMERS = [10, 11, 12]

# DP algorithm adapted from Langmead's notebooks
def _trace(D, x, y):
    ''' Backtrace edit-distance matrix D for strings x and y '''
    i, j = len(x), len(y)
    while i > 0:
        diag, vert, horz = sys.maxsize, sys.maxsize, sys.maxsize
        delt = None
        if i > 0 and j > 0:
            delt = 0 if x[i-1] == y[j-1] else 1
            diag = D[i-1, j-1] + delt
        if i > 0:
            vert = D[i-1, j] + 1
        if j > 0:
            horz = D[i, j-1] + 1
        if diag <= vert and diag <= horz:
            # diagonal was best
            i -= 1; j -= 1
        elif vert <= horz:
            # vertical was best; this is an insertion in x w/r/t y
            i -= 1
        else:
            # horizontal was best
            j -= 1
    # j = offset of the first (leftmost) character of t involved in the
    # alignment
    return j

def _kEditDp(p, t):
    ''' Find the alignment of p to a substring of t with the fewest edits.  
        Return the edit distance and the coordinates of the substring. 
        
        In:
         - p: Pattern, shorter string going along y-axis
         - t: Text, longer string going along x-axis
        
        Out:
         - mn: Minimal number of changes for the pattern to map
         - off: Starting index of mapped pattern
         - mnJ: Ending index of mapped pattern
        '''
    D = numpy.zeros((len(p)+1, len(t)+1), dtype=int)
    # Note: First row gets zeros.  First column initialized as usual.
    D[1:, 0] = range(1, len(p)+1)
    for i in range(1, len(p)+1):
        for j in range(1, len(t)+1):
            delt = 1 if p[i-1] != t[j-1] else 0
            D[i, j] = min(D[i-1, j-1] + delt, D[i-1, j] + 1, D[i, j-1] + 1)
    # Find minimum edit distance in last row
    mnJ, mn = None, len(p) + len(t)
    for j in range(len(t)+1):
        if D[len(p), j] < mn:
            mnJ, mn = j, D[len(p), j]
    # Backtrace; note: stops as soon as it gets to first row
    off = _trace(D, p, t[:mnJ])
    # Return edit distance and t coordinates of aligned substring
    return int(mn), off, mnJ

def partition(word, pieces=2):
    assert len(word) >= pieces
    parts = []
    while pieces > 0:
        part_word_length = (len(word) // pieces)
        parts.append(word[:part_word_length])
        word = word[part_word_length:]
        pieces -= 1
    return parts

# example index
class simpleIndex:
    def __init__(self, text: str):
        self.text = text
        self.kmers = defaultdict(list)

    def setKmers(self, part_length: int):
        text = self.text
        for i in range(len(text) - part_length + 1):
            self.kmers[text[i:i+part_length]].append(i)

    def occurences(self, part: str):
        return self.kmers[part]        

    def query(self, pattern, edist, step):
        hits = []
        for i in range(0, len(pattern)-self.k+1, step):
            for j in self.kmers[pattern[i:i+self.k]]:
                lf = max(0, j-i-edist)
                rt = min(len(self.text), j-i+len(pattern)+edist)
                mn, soff, eoff = _kEditDp(pattern, self.text[lf:rt])
                soff += lf
                eoff += lf
                if mn<=edist:
                    hits.append((mn, soff, eoff))
        hits.sort()
        return hits


    def betterQuery(self, pattern: str):
        text = self.text
        edits = len(pattern) // 10
        pattern_offset = 0
        occurrences = []
        seen = set()
        for part in partition(pattern, edits + 1):
            # print(part)
            for hit in self.occurences(part):
                # print("HIT!")
                lf = max(0, hit - pattern_offset - edits)
                rt = min(len(text), hit - pattern_offset + len(pattern) + edits)
                mn, start_index, end_index = _kEditDp(pattern, text[lf:rt])
                start_index += lf
                end_index += lf
                if mn <= edits and (mn, start_index) not in seen:
                    occurrences.append((mn, start_index, end_index))
                    seen.add((mn, start_index))
            pattern_offset += len(part)
            print("--- %s seconds ---" % (time.time() - start_time))
        return occurrences

# FORMULAS
# len(pattern) // (len(pattern) // 10 + 1)
# if (len(pattern) // (len(pattern) // 10 + 1) != len(text) / (len(pattern) // 10 + 1)):
#     len(pattern) // (len(pattern) // 10 + 1) + 1



seq_rec_list=[seq_record for seq_record in SeqIO.parse(argv[1], "fasta")]
index = simpleIndex(str(seq_rec_list[0].seq))
del seq_rec_list

# Preprocessing including error rate
for i in POSSIBLE_KMERS:
    index.setKmers(i)
print("--- %s seconds ---" % (time.time() - start_time))

fout = open(argv[3], "w")
reads = SeqIO.parse(argv[2], "fasta")
for read in reads:
    hits = index.query(str(read.seq), 100, 1)
    if hits:
        fout.write("{}\t{}\t{}\n".format(read.id, hits[0][1], hits[0][2]))
fout.close()

print("--- %s seconds ---" % (time.time() - start_time))