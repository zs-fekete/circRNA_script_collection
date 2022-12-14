import sys
import os
import argparse
import bamnostic as bs
import pandas as pd
import collections
import numpy as np

sys.path.insert(0, './Script')
import Circstar_readcounts as cr

parser = argparse.ArgumentParser()
#Generic arguments
parser.add_argument("-b", "--bamfile",
                   help = 'merged indexed bamfile of read trios') #all trio .bam files
parser.add_argument("-c", "--circout",
                    help = 'filtered (or at least summed) circSTAR output')

args = parser.parse_args()

#Handling input
bam = bs.AlignmentFile(args.bamfile, 'rb')
circ_data = cr.read_circstar_out(args.circout)
beddata = cr.make_bed(circ_data)

def supporting_reads(bam, bedinfo):
    """Get info about circ-specific supporting reads from bamfile"""
    
    suppreads = []
    
    #Get all reads contained in the circ region
    read_info = [
        (read.query_name, read.reference_start, read.reference_end)
        for read in bam.fetch(contig = bedinfo[0],
                          start = bedinfo[1],
                          end = bedinfo[2])
    ]
    #print(len(read_info))
    
    occ_dic = collections.defaultdict(lambda : [int(0), int(0)])
    for read in read_info:
        occ_dic[read[0]][0] += 1
        if read[1] == bedinfo[1]-1 or read[2] == bedinfo[2]:
            occ_dic[read[0]][1] += 1

    for read in read_info:
        if occ_dic[read[0]][0] == 3 and occ_dic[read[0]][1] >= 2:
            suppreads.append(read[0])
            
    return set(suppreads)

for circid in circ_data.keys():
    bedinfo = beddata[circid]
    reads = list(supporting_reads(bam, bedinfo))

    print(circid, *reads, sep = '\t')
