#! /usr/bin/env python

import pkg_resources
import operator

#PyVCF
import vcf

#BioPython
from Bio import SeqIO

class snpit(object):

    def __init__(self,threshold=10):

        self.threshold=threshold

        # construct the relative path in the package to the library file which contains a list of all the lineages and sub-lineages
        resource_path = '/'.join(('..','lib', 'library'))

        # open a stream object ready for reading
        library_file = pkg_resources.resource_stream("snpit", resource_path)

        self.reference_snps={}

        self.lineages=[]

        for record in library_file:

            lineage_name = record.rstrip().decode('UTF-8')

            self.lineages.append(lineage_name)

            lineage_path='/'.join(('..','lib',lineage_name))

            lineage_file = pkg_resources.resource_stream("snpit", lineage_path)

            self.reference_snps[lineage_name]={}

            for line in lineage_file:

                cols=line.rstrip().decode('UTF-8').split('\t')

                self.reference_snps[lineage_name][int(cols[0])]=cols[1]



    def load_vcf(vcf_file):

        self._reset_lineage_snps()

        vcf_reader = vcf.Reader(open(vcf_file, 'r'))

        for record in vcf_reader:

            for sample in record.samples:
                geno = sample['GT'][0]
                if geno == '.':
                    self._permute(record.POS,'-')
                elif geno != 0:
                    self._permute(record.POS,str(record.ALT[int(geno)-1]))


    def _reset_lineage_snps(self):

        genbank_path = '/'.join(('..','lib', "H37Rv.gbk"))

        # open a stream object ready for reading
        genbank_file = pkg_resources.resource_filename("snpit", genbank_path)

        reference_genome=SeqIO.read(genbank_file,'genbank')

        self.sample_snps={}

        for lineage_name in self.lineages:

            self.sample_snps[lineage_name]={}

            for pos in self.reference_snps[lineage_name]:

                # CAUTION the GenBank File is 1-based, but the lineage files are 0-based
                self.sample_snps[lineage_name][int(pos)]=reference_genome.seq[int(pos)-1]

    def _permute(self,pos,new_base):


        for lineage_name in self.lineages:

            # only proceed if the passed position occurs
            if pos in self.reference_snps[lineage_name].keys():

                # change the
                self.sample_snps[lineage_name][int(pos)]=new_base

    def determine_lineage(self):

        self.percentage={}

        for lineage_name in self.lineages:

            overlap = set(self.reference_snps[lineage_name].items()) & set(self.sample_snps[lineage_name].items())
            shared = float(len(overlap))
            ref = float(len(self.reference_snps[lineage_name]))

            self.percentage[lineage_name]=((shared / ref) * 100)

        results = sorted(self.percentage.items(), key=operator.itemgetter(1),reverse=True)

        # if the first two are above the threshold
        if results[0][1]>self.threshold and results[1][1]>self.threshold:

            # and the first is lineage4, take the second one as it will be a sub-lineage
            if results[0][0]=="lineage4":
                return(results[1])
            # otherwise just go with the highest
            else:
                return(results[0])

        elif results[0][1]>self.threshold:
            return(results[0])

        else:
            return(None)
