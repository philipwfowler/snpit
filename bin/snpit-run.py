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

    tb=snpit(vcf_file=options.input.rstrip(),threshold=10)

    # print("Most likely lineage is ...")
    if tb.percentage is not None:
        print("%s\n%20s %16s %16s %.1f %%" % (options.input,tb.species,tb.lineage,tb.sublineage,tb.percentage))
    else:
        print("%s\n%20s %16s %16s %.1f %%" % (options.input,"Unknown","Unknown",None,0))
