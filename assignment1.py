import mysql.connector
import os
import pybedtools
import pysam

__author__ = 'alexander.bindeus'

class Assignment1:
    
    def __init__(self):
        self.target = "CPT1A"
        self.myGene = self.fetch_gene_coordinates(self, "hg19")

        # hard coded input...do not care
        self.refBamFile = os.path.join(os.getcwd(),
                                       "HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam")

        self.samFile = pysam.AlignmentFile(self.refBamFile, "rb")
        self.reads = list(self.samFile.fetch(self.myGene.chrom, self.myGene.txStart, self.myGene.txEnd))

    @staticmethod
    def fetch_gene_coordinates(self, genome_reference):

        # database has 2 entries labeled as "CPT1A"
        # here i am using the second one

        print("Connecting to UCSC to fetch data of target gene\t", self.target)
        cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password',
                                      db=genome_reference)
        cursor = cnx.cursor()

        query_fields = ["refGene.name2",
                        "refGene.name",
                        "refGene.chrom",
                        "refGene.txStart",
                        "refGene.txEnd",
                        "refGene.strand",
                        "refGene.exonCount",
                        "refGene.exonStarts",
                        "refGene.exonEnds"]

        query = "SELECT DISTINCT %s from refGene" % ",".join(query_fields) + \
                " WHERE refGene.name2=" + '"' + self.target + '"' ""
        cursor.execute(query)

        for row in cursor:
            self.myGene = MyGene(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8])

        cursor.close()
        cnx.close()
        
        print("Done fetching data")

        return self.myGene

    def get_sam_header(self):
        print("\nHeader Info:")
        for key, value in self.samFile.header['HD'].items():
            print('\t', key, value)

    def get_properly_paired_reads_of_gene(self):
        count = 0
        for read in self.reads:
            if read.is_proper_pair:
                count += 1

        print("properly parired end reads", count)

    def get_gene_reads_with_indels(self):
        count = 0
        for read in self.reads:
            if not read.is_unmapped:
                cigarLine = read.cigar
                for (cigarType, cigarLength) in cigarLine:
                    if (cigarType == 1) or (cigarType == 2):  # insertion or deletion
                        count += 1

        print("number of indels:\t", count)

    def calculate_average_coverages(self):
        a = pybedtools.BedTool(self.refBamFile)
        b = a.genome_coverage(bg=True)

        count = 0
        av = 0
        countGene = 0
        avGene = 0

        # slow but seems to work
        for line in b:
            num = float(line[3])
            av += num
            count += 1
            cooBeg = int(line[1])

            # equal sign acutally not needed
            if cooBeg >= self.myGene.txStart:
                if int(line[2]) <= self.myGene.txEnd:
                    avGene += num
                    countGene += 1

        print("average chromosome coverage:\t", av/count)
        print("average gene coverage:\t", avGene/countGene)

    def get_number_mapped_reads(self):
        count = 0
        for read in self.reads:
            if not read.is_unmapped:
                count += 1

        print("unmapped reads", count)

    def get_gene_symbol(self):
        print("symbol for gene:\t", self.myGene.name)

    def get_region_of_gene(self):
        print("region of gene:\n\tchromosome:\t", self.myGene.chrom)
        print("\tcoordinates:\t", self.myGene.txStart, self.myGene.txEnd)

    def get_number_of_exons(self):
        print("number of exons\t", self.myGene.exonCount)

    def print_summary(self):
        print("\nTARGET GENE: ", self.target)
        self.get_sam_header()
        self.get_gene_symbol()
        self.get_properly_paired_reads_of_gene()
        self.get_gene_reads_with_indels()
        self.get_number_of_exons()
        self.get_region_of_gene()
        self.get_number_mapped_reads()
        self.calculate_average_coverages()

        self.samFile.close()

class MyGene:
    def __init__(self, name2, name, chrom, txStart, txEnd, strand, exonCount, exonStarts, exonEnds):
        self.name2 = name2
        self.name = name
        self.chrom = chrom[3:]
        self.txStart = txStart
        self.txEnd = txEnd
        self.strand = strand
        self.exonCount = exonCount
        self.exonStarts = str(exonStarts).lstrip("b'").rstrip(",'").split(",")
        self.exonEnds = str(exonEnds).lstrip("b'").rstrip(",'").split(",")

if __name__ == '__main__':
    print("Assignment 1")
    print(__author__)
    assignment1 = Assignment1()
    assignment1.print_summary()
