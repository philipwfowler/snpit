#! /usr/bin/env python

from snpit import snpit

import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf",required=True,help="the path to a VCF file ")
    parser.add_argument("--genbank",default="H37Rv.gbk",help="the path to the genbank file of the H37Rv M. tuberculosis reference collection genome")
    options = parser.parse_args()

    # create an instance (this loads all the lineages)
    tb_lineage_collection=snpit(threshold=10)

    # load the VCF file into the SNPIT instance (this automatically resets the list of SNPs for each lineage)
    tb_lineage_collection.load_vcf(options.vcf.rstrip())

    # determine the most likely lineage/sub-lineage using snpit
    (lineage,percentage)=tb_lineage_collection.determine_lineage()


    print("Most likely lineage is ...")
    if percentage is not None:
        print("%16s %.1f %%" % (lineage,percentage))
    else:
        print("%16s %.1f %%" % ("None identified",0))

    print("\n")
    print("**************************")
    print("List of lineages identified in descending order")

    for (lineage,percentage) in tb_lineage_collection.results:
        print("%16s %.1f %%" % (lineage,percentage))
