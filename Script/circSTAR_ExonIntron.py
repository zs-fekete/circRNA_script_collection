import sys
import argparse
import os
import collections
import bamnostic as bs
import pandas as pd
import numpy as np


import Circstar_readcounts as cr
import Circstar_bamtrios as cb


parser = argparse.ArgumentParser()
#Generic arguments
parser.add_argument("-b", "--bamfile",
                   help = 'merged indexed bamfile of read trios') #all trio .bam files
parser.add_argument("-c", "--circout",
                    help = 'filtered (or at least summed) circSTAR output')
parser.add_argument("-g", "--gtf_file", help = 'gtf file to use for annotation')
parser.add_argument("-o", "--out_prefix", help = "Output file prefix") #base values are set at the relevant function

#Optional arguments
parser.add_argument("--ei_limit", type = float, default = 0.2,
                   help = 'intron/exon coverage to accept an EIciRNA')

args = parser.parse_args()


gtf_data = cb.read_gtf(args.gtf_file)
circ_data = cr.read_circstar_out(args.circout)
bam = bs.AlignmentFile(args.bamfile, 'rb')

cb.print_eicirna_summary(circ_data, bam, gtf_data, args.ei_limit, base = args.out_prefix)