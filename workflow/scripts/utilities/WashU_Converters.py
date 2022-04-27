import sys
import json
import argparse
from itertools import chain

# Conversions available are:
# 1) bed intervals to bed intervals (--convert-from bp --convert-to vcf)
# 2) single basepair file to variant file (--convert-from bp --convert-to vcf)
# 3) bedpe intervals to longrange (--convert-from bedpe --convert-to longrange)

# # Making a command line interface
parser = argparse.ArgumentParser()
parser.add_argument('--convert-from', type=str, choices=['bed', 'bedpe', 'bp'], required=True)
parser.add_argument('--convert-to', type=str, choices=['bed', 'longrange', 'vcf'], required=True)
parser.add_argument('--input-chr-format', action='store_true', help='Specify the presence of chr in the chromosome fields.')
parser.add_argument('--output-chr-format', action='store_true', help='Specify the output of chr in the chromosome fields.')

# arguments for bed 
parser.add_argument('--chr', type=int, default=1)
parser.add_argument('--start', type=int, default=2)
parser.add_argument('--end', type=int, default=3)

# arguments for bedpe 
parser.add_argument('--chrB', type=int, default=4)
parser.add_argument('--startB', type=int, default=5)
parser.add_argument('--endB', type=int, default=6)
parser.add_argument('--score', type=int, default=7)

# arguments for vcf
parser.add_argument('--pos', type=int, default=-1)
parser.add_argument('--id', type=int, default=None)
parser.add_argument('--ref', type=int, default=None)
parser.add_argument('--alt', type=int, default=None)
parser.add_argument('--qual', type=int, default=None)
parser.add_argument('--filter', type=int, default=None)
parser.add_argument('--info', type=int, default=None)
parser.add_argument('--format', type=int, default=None)
parser.add_argument('--samples', type=int, default=None)

params = parser.parse_args()

# set the clean_chr function for parsing the chr fields 
if params.input_chr_format == True and params.output_chr_format == True:
    clean_chr = lambda x: x

elif params.input_chr_format == True and params.output_chr_format == False:
    clean_chr = lambda x: x.replace('chr', '')

elif params.input_chr_format == False and params.output_chr_format == True:
    clean_chr = lambda x: 'chr{}'.format(x)

elif params.input_chr_format == False and params.output_chr_format == False:
    clean_chr = lambda x: x


#########################################################################################
# Process a bp into a vcf file
#########################################################################################

if params.convert_from == 'bp' and params.convert_to == 'vcf':

    # convert to zero-based indexing
    params.chr -= 1
    params.pos -= 1

    # write out the header
    vcf_header = "##fileformat=VCFv4.3\n"
    vcf_header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    print(vcf_header)

    # parse the file 
    for i, line in enumerate(sys.stdin):

        line_info = line.strip().split('\t')
        chrom = line_info[params.chr].replace('chr', '')
        pos = line_info[params.pos]
        idd = line_info[params.id] if params.id != None else '.'
        ref = line_info[params.ref] if params.ref != None else '.'
        alt = line_info[params.alt] if params.alt != None else '.'
        qual = line_info[params.qual] if params.qual != None else '.'
        flt = line_info[params.filter] if params.filter != None else 'PASS'
        info = life_info[params.info] if params.info != None else '.'
        #form = line_info[params.format] if params.format != None else '.'
        #samples = line_info[params.samples] if params.samples != None else '.'

        new_info = [chrom, pos, idd, ref, alt, qual, flt, info] #, form, samples]
        new_info = '\t'.join(new_info)
        print(new_info)
    

#########################################################################################
# Process a bedpe into a longrange file
#########################################################################################

if params.convert_from == 'bedpe' and params.convert_to == 'longrange':

    # convert to zero-based indexing
    params.chr -= 1
    params.start -= 1
    params.end -= 1
    params.chrB -= 1
    params.startB -= 1
    params.endB -= 1
    params.score -= 1


    # parse the file 
    for i, line in enumerate(sys.stdin):

        # extracting the main fields
        line_info = line.strip().split('\t')

        chrom = clean_chr(line_info[params.chr])
        start = line_info[params.start]
        end = line_info[params.end]
        chromB = clean_chr(line_info[params.chrB])
        startB = line_info[params.startB]
        endB = line_info[params.endB]
        score = line_info[params.score]

        # print out as tab delimited data
        second_anchor = '{}:{}-{},{}'.format(chromB, startB, endB, score)
        new_info = [chrom, start, end, second_anchor]
        new_info = '\t'.join(new_info)
        print(new_info)
