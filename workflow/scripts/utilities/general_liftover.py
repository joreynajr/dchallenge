import os 
import argparse
from numpy import nan
from biolib import liftover as liftover_utility
from liftover import get_lifter
import codecs
def unescaped_str(arg_str):
        return codecs.decode(str(arg_str), 'unicode_escape')

lift_converter = get_lifter('hg38', 'hg19')

# commandline interface
parser = argparse.ArgumentParser()
parser.add_argument('-i', type=str, help='Path to the input file.')
parser.add_argument('-o', type=str, help='Path to the output file.')
parser.add_argument('--chr-col', type=int, help='One-based column index for the chromosome.')
parser.add_argument('--pos-col', type=int, help='One-based column index for the position.')
parser.add_argument('--header', action='store_true', help='Boolean indicating the presense of a header.')
parser.add_argument('--sep', type=unescaped_str, help='Separator used.')
params = parser.parse_args()

#fn = 'results/main/coloc/Results/eQTL_Catalogue/T1D_34012112_Gaulton/BLUEPRINT/T-cell/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed'
#outfn = 'test.tsv'
chr_col = params.chr_col - 1
pos_col = params.pos_col - 1
header = True

#print("params: ", params)

with open(params.i) as fr, open(params.o, 'w') as fw:
    
    # write the header
    if params.header:
        header_str = next(fr).strip()
        write_str = params.sep.join([header_str, 'old_pos', 'was_converted'])
        fw.write(write_str + '\n')
        
    # parse the rest of the file
    for line in fr:
        
        # get chrom and pos
        info = line.strip().split(params.sep)
        #print('info: ', info)
        chrom = info[chr_col]
        pos = int(info[pos_col])
        
        # get the new chrom and pos
        new_info = info.copy()
        new_coords = lift_converter[chrom][pos]
                 
        if len(new_coords) == 0:
            new_info[chr_col] = nan
            new_info[pos_col] = nan
            new_info.extend([pos, 0])
            
        elif len(new_coords) == 1:
            new_info[chr_col] = new_coords[0][0]
            new_info[pos_col] = new_coords[0][1]
            new_info.extend([pos, 1])

        else:
            raise Exception('Check this case: {}'.format(new_coords))
  
        new_info = [str(x) for x in new_info]
        new_info_str = params.sep.join(new_info) + '\n'
        fw.write(new_info_str)
