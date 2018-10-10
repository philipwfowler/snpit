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
    parser.add_argument("--input",required=True,help="the path to the VCF or FASTA file to read and classify (either can be bzip2ed/gzipped)")
    options = parser.parse_args()

    # now try opening the specified input file
    try:
        TEST = open(options.input,'r')
    except IOError:
        print("input file "+options.input+" does not exist!")

    # create an instance (this loads all the lineages)
    tb=snpit(threshold=10,input_file=options.input)

    if tb.percentage is not None:
        print("%s\n%16s %16s %16s %.1f %%" % (options.input,tb.species,tb.lineage,tb.sublineage,tb.percentage))
    else:
        print("%s\n%16s" % (options.input,"none identified"))
