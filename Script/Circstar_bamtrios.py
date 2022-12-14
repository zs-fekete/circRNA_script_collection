import sys
import os
import bamnostic as bs
import pandas as pd
import collections
import numpy as np

import Circstar_readcounts as cr

#This is explicitly for stuff supporting the whole circRNA
#does NOT work for reads partially overlapping a feature
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
    print(len(read_info))
    
    occ_dic = collections.defaultdict(lambda : [int(0), int(0)])
    for read in read_info:
        occ_dic[read[0]][0] += 1
        if read[1] == bedinfo[1]-1 or read[2] == bedinfo[2]:
            occ_dic[read[0]][1] += 1

    for read in read_info:
        if occ_dic[read[0]][0] == 3 and occ_dic[read[0]][1] >= 2:
            suppreads.append(read)
            
    return suppreads

#For what the previous won't handle
def feature_overlap_reads(bam, bedinfo):
    """Get info about reads overlapping a feature within circ"""
    suppreads = []
    
    #Get all reads contained in the circ region
    read_info = [
        (read.query_name, read.reference_start, read.reference_end)
        for read in bam.fetch(contig = bedinfo[0],
                          start = bedinfo[1],
                          end = bedinfo[2])
    ]
            
    return read_info



#Read gtf and get features from it
def read_gtf(filename):
    
    gtf_df = pd.read_csv(filename,
           sep = '\t',
           names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'],
           comment = '#')
    
    return(gtf_df)

#the get_exons and get_introns functions return the location in bed style, tuple type
def get_exons(bedinfo, gtf_df):
    """Get exons that are annotated in the gtf"""
    exons = []
    
    exon_overlap = gtf_df[
    (gtf_df['seqname'] == bedinfo[0]) &
    (gtf_df['feature'] == 'exon') &
    ((gtf_df['end'].between(bedinfo[1], bedinfo[2], inclusive = 'both')) |
    (gtf_df['start'].between(bedinfo[1], bedinfo[2], inclusive = 'both'))|
    ((gtf_df['start'] <= bedinfo[1]) & (gtf_df['end'] >= bedinfo[2])))
    ]
    
    for index, exon in exon_overlap.iterrows():
        exons.append((exon['seqname'], exon['start'], exon['end']))
        
    return exons

def get_introns(bedinfo, exons):
    """Get introns aka everything not overlapping with an exon"""
    
    all_pos = range(bedinfo[1], bedinfo[2])
    exons = list(get_coverage(exons).keys())
    
    intronic = [i for i in all_pos if i not in exons]
    
    if len(intronic) > 0:
        return get_intron_ranges(intronic, bedinfo[0])
    else:
        return [] #this prevents a NoneType error later
        

def get_intron_ranges(intronic, contig):
    """Get all intronic positions"""
    
    starts = [ intronic[i] for i in range(len(intronic)-1) 
              if intronic[i+1] != intronic[i]+1 ]
    starts.insert(0, intronic[0])
    ends = [ intronic[i+1] for i in range(len(intronic)-1) 
            if intronic[i+1] != intronic[i]+1 ]
    ends.append(intronic[-1])
    
    introns = [(contig, starts[i], ends[i]) for i in range(0, len(starts))]
    
    return introns


#Coverages and defining EIciRNAs
def get_coverage(suppreads):
    """Getting per coordinate coverage for the whole circ interval"""
    
    #positions = range(bedinfo[1]-1, bedinfo[2])
    readcov = [list(range(read[1], read[2])) for read in suppreads]
    coverages = collections.Counter([point for read in readcov for point in read])
    
    if len(coverages) == 0:
        return {"mock" : 0}
    else:
        return coverages


def feature_type_sum(features, bam):
    feature_covs = [ np.mean(list(get_coverage(feature_overlap_reads(bam, f)).values()))
        for f in features ]
        
    return feature_covs
    

def is_eicirna(exons, introns, limit, bam):
    if len(exons) == 0:
        intron_cov = str(np.mean(feature_type_sum(introns, bam)))
        return ['0', intron_cov, 'intronic']
    if len(introns) == 0:
        exon_cov = str(np.mean(feature_type_sum(exons, bam)))
        return [exon_cov, '0', 'exonic']
    else:
        exon_cov = str(np.mean(feature_type_sum(exons, bam)))
        intron_cov = str(np.mean(feature_type_sum(introns, bam)))
        if float(exon_cov) == 0:
            return [exon_cov, intron_cov, 'intronic']
        elif float(intron_cov) / float(exon_cov) > limit:
            return [exon_cov, intron_cov, 'EIciRNA']
        else:
            return [exon_cov, intron_cov, 'exonic/unlikely EIciRNA']
    
    
    
#Gene types and names
#(Might move to other script later)
def get_genes(bedinfo, gtf_df):
    """Get everything about genes overlapping with circRNA"""
    genes = []
    all_pos = set(range(bedinfo[1], bedinfo[2]))
    
    gene_overlap = gtf_df[
    (gtf_df['seqname'] == bedinfo[0]) &
    (gtf_df['feature'] == 'gene')
    ]
    
    for index, gene in gene_overlap.iterrows():
        generange = set(range(gene['start'], gene['end']))
        
        #keep it if found in both sets (= there is an overlap)
        if (all_pos & generange):
            genes.append(gene)

    return genes


def get_gene_type(gene):
    
    genetype = gene[8].split('gene_biotype "')[1].split('"')[0]
    return genetype

def all_gene_types(genes):
    if len(genes) > 0:
        types = [get_gene_type(gene) for gene in genes]
        return types
    else:
        return []
    
def get_gene_id(gene):
    
    genename = gene[8].split('gene_id "')[1].split('"')[0]
    return genename

def all_gene_ids(genes):
    if len(genes) > 0:
        names = [get_gene_id(gene) for gene in genes]
        return names
    else:
        return []
    
    
#Output stuff
def print_eicirna_summary(circ_data, bam, gtf_data, ei_limit, base = 'exon_intron'):
    """Write EI-annotation results to a tsv"""
    
    output = open(f'{base}.tsv', "w")
    output.write(f'circID\tgene_types\tgene_ids\texon_coverage\tintron_coverage\tcirc_type\n')
    
    for circid in circ_data.keys():
        beddata = cr.make_bed(circ_data)[circid]
        genes = get_genes(beddata, gtf_data)
        exons = get_exons(beddata, gtf_data)
        introns = get_introns(beddata, exons)
        circtype = is_eicirna(exons, introns, ei_limit, bam)
        
        output.write(f'{circid}\t')
        if len(genes) > 0:
            output.write(';'.join(all_gene_types(genes))+'\t')
            output.write(';'.join(all_gene_ids(genes))+'\t')
        else:
            output.write(f'No annotation\tNo annotation\t')
        output.write('\t'.join(circtype)+'\n')

    output.close