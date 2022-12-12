import sys
import argparse
import os
import subprocess
import Circstar_readcounts as cr

parser = argparse.ArgumentParser()
#Generic arguments
parser.add_argument("-i", "--infiles", nargs = '*') #all circstar output .txts
parser.add_argument("-o", "--out_prefix", help = "Output file prefix") #base values are set at the relevant function

#Filtering arguments
parser.add_argument("--max_genomic_len", type = int, default = 0,
                    help = "maximum distance between BSJ junction sites for acceptable circRNA. \
                    Recommended: longest known gene in your organism \
                    Consider: plants may produce 1 circRNA using multiple genes")
parser.add_argument("--min_readpair", type = int, default = 15,
                   help = "Minimum amount of total read pairs supporting the circRNA, across samples")
parser.add_argument("--min_occurrence", type = int, default = 3,
                   help = "Minimum amount of samples where the circRNA is found")
parser.add_argument("--min_readpair_per_sample", type = int, default = 1,
                   help = "Minimum amount of supporting read pairs per sample")

parser.add_argument("--print_unfiltered", action = 'store_true',
                   help = "Print unfiltered, combined data if true")

args = parser.parse_args()

#get inputs
circouts = list(args.infiles)
cirdata = [cr.read_circstar_out(f) for f in circouts]

print("Inferred sample names:")
for i in cr.get_sample_names(circouts):
    print(i)
print()

print(f"Filtering parameters: \n\
      max genomic length: {args.max_genomic_len}, \n\
      minimum read pairs per circRNA: {args.min_readpair} , \n\
      found in at least {args.min_occurrence} samples, \n\
      with at least {args.min_readpair_per_sample} reads per sample")
if args.print_unfiltered:
    print("Also printing unfiltered data")
    
#Run basic summary and filtering 
readcount_summary = cr.summarize_circstar_outs(cirdata)
if args.print_unfiltered:
    cr.print_circstar_summary(readcount_summary,
                              circouts = circouts,
                              base = 'unfiltered_summary')
    cr.output_bed(cr.make_bed(readcount_summary),
                              base = 'unfiltered')
    print(f"Unfiltered circ read data printed to \n\t unfiltered_summary.tsv \n\t unfiltered.bed")
    

filtered_data = cr.filter_circread_amount(
    readcount_summary,
    args.min_occurrence,
    args.min_readpair,
    args.min_readpair_per_sample
)

if args.max_genomic_len > 0:
    filtered_data = cr.filter_circ_length(filtered_data,
                                          args.max_genomic_len)

#Printing output
if args.out_prefix:
    cr.print_circstar_summary(filtered_data,
                              circouts = circouts,
                              base = args.out_prefix)
    cr.output_bed(cr.make_bed(filtered_data),
                  base = args.out_prefix)
    print(f"Output printed to: \n\t {args.out_prefix}.tsv \n\t {args.out_prefix}.bed")
else:
    cr.print_circstar_summary(filtered_data,
                             circouts = circouts)
    cr.output_bed(cr.make_bed(filtered_data))
    print(f"Output printed to: \n\t circ_readpair_summary.tsv \n\t circ_summary.bed")


          