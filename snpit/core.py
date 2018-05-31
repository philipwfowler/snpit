#! /usr/bin/env python

import pkg_resources
import operator

import vcf

#BioPython
from Bio import SeqIO

class snpit(object):

    """
    The snpit class is designed to take a VCF file and return the most likely lineage based on Sam Lipworth's SNP-IT.

    The methods have been separated so it can be incorporated into single Python scripts that processes multiple VCF files.
    """

    def __init__(self,threshold=10):

        """
        Args:
            threshold: The percentage of snps above which a sample is considered to belong to a lineage.
        """

        # set the threshold as a class attribute
        self.threshold=threshold

        # construct the relative path in the package to the library file which contains a list of all the lineages and sub-lineages
        resource_path = '/'.join(('..','lib', 'library'))

        # open a stream object ready for reading
        library_file = pkg_resources.resource_stream("snpit", resource_path)

        self.reference_snps={}

        self.lineages=[]

        # read the library file line-by-line
        for record in library_file:

            # remove the carriage return and decode from binary
            lineage_name = record.rstrip().decode('UTF-8')

            # remember the lineage name in a list
            self.lineages.append(lineage_name)

            # now we know the name construct the relative path to this lineage file
            lineage_path='/'.join(('..','lib',lineage_name))

            # open a stream object to that file ready for reading
            lineage_file = pkg_resources.resource_stream("snpit", lineage_path)

            # initialise the dictionary for this lineage
            self.reference_snps[lineage_name]={}

            # read the lineage file, line-by-line
            for line in lineage_file:

                # remove the carriage return, decode from binary, and split on tabs
                cols=line.rstrip().decode('UTF-8').split('\t')

                # remember the base in the dictionary using the genome position as the key
                self.reference_snps[lineage_name][int(cols[0])]=cols[1]


    def load_vcf(self,vcf_file):
        """
        Loads the vcf file and then, for each lineage, identify the base at each of the identifying positions in the genome.

        Args:
            vcf_file: Path to the VCF file to be read
        """

        # setup the dictionaries of expected SNPs for each lineage
        self._reset_lineage_snps()

        # open the VCF file for reading
        vcf_reader = vcf.Reader(open(vcf_file, 'r'))

        # read the VCF file line-by-line
        for record in vcf_reader:

            # consider each lineage in turn
            for lineage_name in self.lineages:

                # only proceed if the genome position occurs in the list of identifiable positions
                if record.POS in self.reference_snps[lineage_name].keys():

                    # parse the record
                    for sample in record.samples:
                        geno = sample['GT'][0]

                        # if there is a null call, record a hyphen which won't match, regardless of the reference
                        if geno == '.':
                            self.sample_snps[lineage_name][int(pos)]="-"
                            self._permute(record.POS,'-')

                        # otherwise replace the H37Rv base with the actual in the VCF file
                        elif geno != 0:

                            # CAUTION the GenBank File is 1-based, but the lineage files are 0-based
                            self.sample_snps[lineage_name][int(pos)]=record.ALT[int(geno)-1]

    def _reset_lineage_snps(self):
        """
        For each lineage creates a dictionary of the positions and expected nucleotides for TB that
        define that lineage.

        This is required because the VCF files only list changes relative to H37Rv.
        Hence these dictionaries are then changed when mutations at these positions are encountered.
        """

        # make the relative path to the H37Rv TB reference GenBank file
        genbank_path = '/'.join(('..','lib', "H37Rv.gbk"))

        # open a stream object ready for reading
        genbank_file = pkg_resources.resource_filename("snpit", genbank_path)

        # read the reference genome using BioPython
        reference_genome=SeqIO.read(genbank_file,'genbank')

        self.sample_snps={}

        #  iterate through the lineages
        for lineage_name in self.lineages:

            self.sample_snps[lineage_name]={}

            # iterate over the positions in the reference set of snps for that lineage
            for pos in self.reference_snps[lineage_name]:

                # CAUTION the GenBank File is 1-based, but the lineage files are 0-based
                # Remember the nucleotide at the defining position
                self.sample_snps[lineage_name][int(pos)]=reference_genome.seq[int(pos)-1]



    def determine_lineage(self):
        """
        Having read the VCF file, for each lineage, calculate the percentage of SNP present in the sample.
        Note that this means the percentages will not add up to 100%.

        Returns:
            tuple of (lineage,percentage)
        """

        self.percentage={}

        # consider lineage-by-lineage
        for lineage_name in self.lineages:

            # using sets, calculate which SNPs at the defining positions are the same as the reference for this lineage
            overlap = set(self.reference_snps[lineage_name].items()) & set(self.sample_snps[lineage_name].items())

            # calculate how many there are
            shared = float(len(overlap))

            # calculate how many SNPs define this lineage
            ref = float(len(self.reference_snps[lineage_name]))

            # thereby calculate the percentage of SNPs in this sample that match the lineage
            self.percentage[lineage_name]=((shared / ref) * 100)

        # create an ordered list of tuples of (lineage,percentage) in descending order
        results = sorted(self.percentage.items(), key=operator.itemgetter(1),reverse=True)

        # if the first two entries are above the threshold
        if results[0][1]>self.threshold and results[1][1]>self.threshold:

            # and the first is lineage4, take the second one as it will be a sub-lineage
            if results[0][0]=="lineage4":
                return(results[1])

            # otherwise just return the highest
            else:
                return(results[0])

        # otherwise if the first element is above the threshold
        elif results[0][1]>self.threshold:

            # return the tuple
            return(results[0])

        # finally, no strain must be above the threshold percentage
        else:
            return((None,None))
