#! /usr/bin/env python

#programme to report the lineage/subspecies of a TB sample alligned to NC000962
#useage: python SNP-IT.py [guuid] [name of outfile]
#or for lots of samples: cat [list of samples] | parallel -j[no. of threads] python fasta_typer9.py {} {}.out
#Output1: is to file containing absolute and % hits for all subspecies
#Output2: only the top call, is to standard out - redirect to file if you want eg 1>calls.log
import vcf

from snpit import snpit

import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf",required=False,help="the path to the VCF file to read and classify (can be zip/gzipped)",dest='vcf')
    parser.add_argument("--fasta",required=False,help="the path to the fasta file (gzipped)",dest='fasta')
    options = parser.parse_args()

    fasta = options.fasta
    vcf= options.vcf

    if fasta:
        input = fasta
    else:
        input =vcf
    print(input)
    # now try opening the specified input file
    try:
        TEST = open(input,'r')
    except IOError:
        print("input file "+input+" does not exist!")


    # create an instance (this loads all the lineages)
    tb_lineage_collection=snpit(threshold=10)

    # load the VCF file into the SNPIT instance (this automatically resets the list of SNPs for each lineage)
    if vcf:
    	tb_lineage_collection.load_vcf(input.rstrip())
    elif fasta:
    	tb_lineage_collection.load_fasta(input.rstrip())
    else:
        sys.stdout.write('Please provide input files')
    # determine the most likely lineage/sub-lineage using snpit
    (lineage,percentage)=tb_lineage_collection.determine_lineage()

    # print("Most likely lineage is ...")
    if percentage is not None:
        	print("%40s %16s %.1f %%" % (input,lineage,percentage))
    else:
        print("%40s %16s %.1f %%" % (input,"none identified",0))
