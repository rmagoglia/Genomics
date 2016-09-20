"""switch_chrom_names.py

Takes a text file and replaces chromosome names with the appropriate names
from a different naming convention. Supported formats:

UCSC: chr1, chr2, ...
RefSeq: NC_000001.11, NC_000002.12, ...
Ensembl: 1, 2, ...

Note that converting between UCSC and Ensembl formats can be done with a 
simple awk command; this script currrently only converts between refseq 
and the other two formats. 

"""

import argparse

def get_chrom_dict(chrom_file, input_format, output_format):
    """Reads chromosomes into a dictionary."""

    chromDict = {}

    with open(chrom_file, 'r') as chroms:
        for line in chroms:
            refseq, ens = line.strip().split('\t')
            ucsc = "chr" + ens
            if input_format == output_format:
                print("Please select the desired conversion formats.")
                break

            if input_format == "refseq":
                if output_format == "ensembl":
                    chromDict[refseq] = ens
                elif output_format == "ucsc":
                    chromDict[refseq] = ucsc

            elif output_format == "refseq":
                if input_format == "ensembl":
                    chromDict[ens] = refseq
                elif input_format == "ucsc":
                    chromDict[ucsc] = refseq

    return chromDict


def convert_chroms(chrom_dict, input_file, output_file):
    """Reads input file and writes output file with name conversions."""

    converted = 0
    unconverted = 0

    out = open(output_file, 'w')

    with open(input_file, 'r') as inFile:
        for line in inFile:
            chrom = line.strip().split('\t')[0]
            fields = line.strip().split('\t')[1:]
            if chrom in chrom_dict:
                newChrom = chrom_dict[chrom]

                outlist = [newChrom] + fields
                output = '\t'.join(outlist)
                out.write(output + '\n')
                converted += 1
            else:
                unconverted +=1

    print("Total entires: " + str(converted + unconverted))
    print("  Converted: " + str(converted))
    print("  Unconverted: " + str(unconverted))

    out.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("input_file", help="File which needs name conversion, "
        "chromosome names must be the first column, must be tab delimited.")
    parser.add_argument("output_file", help="File with converted names.")
    parser.add_argument("chromosomes_file", help="File containing refseq and "
        "ensembl chromosome names. E.g. NC_000001.11     1")
    parser.add_argument("-a", "--input_format", default="refseq", 
        dest="input_format", help="Chromosome format of input file. Default: refseq")
    parser.add_argument("-b", "--output_format", default="ensembl", 
        dest="output_format", help="Chromosome format of output file. Default: ensembl")

    options = parser.parse_args()

    CHROM_DICT = get_chrom_dict(options.chromosomes_file, options.input_format.lower(),
        options.output_format.lower())

    convert_chroms(CHROM_DICT, options.input_file, options.output_file)
