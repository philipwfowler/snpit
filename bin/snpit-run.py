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
    parser.add_argument("--input",required=True,help="the path to the VCF file to read and classify (can be zip/gzipped)")
    options = parser.parse_args()

    # now try opening the specified input file
    try:
        TEST = open(options.input,'r')
    except IOError:
        print("input file "+options.input+" does not exist!")


    # create an instance (this loads all the lineages)
    tb_lineage_collection=snpit(threshold=10)

    # load the VCF file into the SNPIT instance (this automatically resets the list of SNPs for each lineage)
    tb_lineage_collection.load_vcf(options.input.rstrip())

    # determine the most likely lineage/sub-lineage using snpit
    (lineage,percentage)=tb_lineage_collection.determine_lineage()

    # print("Most likely lineage is ...")
    if percentage is not None:
        print("%16s %.1f %%" % (lineage,percentage))
    else:
        print("%16s %.1f %%" % ("None identified",0))

    print(tb_lineage_collection.results)
