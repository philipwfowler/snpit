#! /usr/bin/env python

import pkg_resources, csv, codecs
import operator

# PyVCF
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
        resource_path = '/'.join(('..','lib', 'library.csv'))

        # open a stream object ready for reading
        library_file = pkg_resources.resource_stream("snpit", resource_path)

        utf8_reader = codecs.getreader("utf-8")


        self.reference_snps={}

        self.lineages={}

        # with open(library_file) as CSVFILE:

        reader = csv.DictReader(utf8_reader(library_file))

        # read the library file line-by-line
        for record in reader:

            # remove the carriage return and decode from binary
            lineage_name = record['id']

            # remember the lineage meta data in a dictionary
            self.lineages[lineage_name]={'species':record['species'],'lineage':record['lineage'],'sublineage':record['sublineage']}

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

                            self.sample_snps[lineage_name][int(record.POS)]="-"

                        # otherwise replace the H37Rv base with the actual base from the VCF file
                        elif geno != 0:
                            self.sample_snps[lineage_name][int(record.POS)]=record.ALT[int(geno)-1]

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

        # consider lineage-by-lineage
        for lineage_name in self.lineages:

            reference_set=[]

            shared=0
            ref=0

            for i,j in enumerate(self.reference_snps[lineage_name]):

                if self.reference_snps[lineage_name][j] == self.sample_snps[lineage_name][j]:
                    shared+=1
                ref+=1

            # thereby calculate the percentage of SNPs in this sample that match the lineage
            self.percentage[lineage_name]=((shared / ref) * 100)

        # create an ordered list of tuples of (lineage,percentage) in descending order
        self.results = sorted(self.percentage.items(), key=operator.itemgetter(1),reverse=True)

        identified_lineage_name=self.results[0][0]
        identified_lineage_percentage=self.results[0][1]

        # if the top lineage is above the specified threshold, return the classification
        if identified_lineage_percentage>self.threshold:
            return(self.lineages[identified_lineage_name]['species'],self.lineages[identified_lineage_name]['lineage'],self.lineages[identified_lineage_name]['sublineage'],identified_lineage_percentage)

        # finally, no strain must be above the threshold percentage so return Nones as "Don't know"
        else:
            return(None,None,None,None)
