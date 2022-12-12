import collections

#Reading in circSTAR outputs

def read_circstar_out(circstar_out):
    """Read in output files of circstar.pl"""
    
    with open(circstar_out) as f:
        cirdic = ({line.split()[0] : line.split()[1:] for line in f.readlines()[1:]})
        
    return cirdic

def get_sample_names(circouts):    
    """Create meaningful sample names from filename"""
    
    #circouts should be the input filenames
    samples = [i.split("/")[-1].split(".")[0] for i in circouts]
    
    return samples


#Summarizing circRNA info across samples

def summarize_circstar_outs(cirdic_list):
    """Make a dictionary of all circs and read counts"""
    
    #Make a set of all circRNA IDs
    circid_list = [list(cirdic.keys()) for cirdic in cirdic_list]
    circid_set = set([item for sublist in circid_list for item in sublist])

    circ_counts = collections.defaultdict(lambda: [])
    #Check the occurrence of each circRNA ID in each sample and add read counts to list
    for circid in circid_set:
        for dic in cirdic_list:
            if circid in dic.keys():
                #append only total read counts, not all info
                circ_counts[circid].append(dic[circid][0])
            else:
                circ_counts[circid].append(0)

    return circ_counts


def make_bed(circ_dict):
    """Get location info from any circ dict"""
    
    circbed_data = collections.defaultdict(lambda : ())
    
    for cir in list(circ_dict.keys()):        
        chrom = cir.split(":")[0]
        try:
            pos1 = int(cir.split(":")[1].split('-')[0])
            pos2 = int(cir.split("-")[1])
        except IndexError:
            if cir == "Region":
                pass
            else:
                print(f"IndexError in make_bed, the culprit was: {cir}")
        
        feature_len = abs(pos1 - pos2)
        
        if pos1 < pos2:
            circbed_data[cir] = (chrom, pos1, pos2, feature_len)
        else:
            circbed_data[cir] = (chrom, pos2, pos1, feature_len)
        
    return circbed_data


#Filtering steps

def filter_circread_amount(circ_counts, min_n = 3, min_val = 15, per_sample = 5):
    """Filter read for read count and occurrence"""
    
    filtered_dic = collections.defaultdict(lambda: [])
    
    for circid, counts in circ_counts.items():
        n = 0
        a = 0
        sumtotal = 0
        for c in counts: #make these with counters?
            try:
                if int(c) >= per_sample: #passes per sample minimum
                    a += 1
                if int(c) > 0: #it's there
                    n += 1
                    sumtotal = sumtotal + int(c)
            except ValueError:
                print(f"Warning: {c} is not an int, skipping")
                    
        #see if it passes all criteria
        if n >= min_n and a >= min_n and sumtotal >= min_val:
            filtered_dic[circid] = counts
    
    return filtered_dic

def filter_circ_length(any_circ_dict, max_len = 0):
    """Filter for genomic length"""

    positions = make_bed(any_circ_dict)
    length_passed = collections.defaultdict(lambda : 0)
    
    for circ in positions.keys():
        if max_len == 0 or int(positions[circ][3]) <= max_len:
            length_passed[circ] = any_circ_dict[circ]
                              
    return length_passed
    

#Printing summary file(s)

def output_bed(beddata, base = 'circ_summary'):
    """Write an output file in bed format"""
    
    output = open(f'{base}.bed', 'w')
    for circ, info in beddata.items():
        for i in info[:-1]:
            output.write(str(i)+'\t')
        output.write(f"{circ}\n")
    output.close

    
def print_circstar_summary(readdata, circouts, base = "circ_readpair_summary"):
    """Write a tsv summary with read counts from any circ dict"""
    
    #Print header
    output = open(f'{base}.tsv', 'w')
    output.write(f'Region\t')
    for sample in get_sample_names(circouts):
        output.write(sample+'\t')
    output.write(f'Total_count\tOccurrence\n')
    
    #Print body
    for circ, info in readdata.items():
        output.write(f"{circ}\t")
        for i in info:
            output.write(str(i)+'\t')
        output.write(f"{sum(map(int, info))}\t")
        output.write(str(sum(int(i) > 0 for i in info))+'\n')
    output.close
