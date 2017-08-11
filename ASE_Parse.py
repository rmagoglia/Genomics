"""ASE_Parse.py 

Takes ASEr output from hybrid and parental data sets
and integrates the results to generate lists of
filtered genes for GO enrichment.

Goals:
 - GO list of hybrid genes for which parental data 
 agrees in directionality
 - GO list of the corresponding parental data

Ouput files:
 - Hybrid GO
 - Parental GO

 Input files:
 - Hybrid files mapped on parent 1 genome
 - Hybrid files mapped on parent 2 genome
 - Parental files mapped on corresponding genomes
 """

from argparse import ArgumentParser 
from scipy.stats import binom_test
from math import log2


def data_to_dict(dataFile):
    '''Reads ASEr output data and returns it as a dict with genes
    as keys and ref / alt expression as values.'''

    data = {}

    with open(dataFile, 'r') as dataFile:
        librarySize = dataFile.readline().split()[-1]
        data['librarySize'] = librarySize

        dataFile.readline() # header line

        for line in dataFile:
            line = line.strip().split()
            gene = line[0]
            info = line[2:4] # ref and alt counts
            data[gene] = info

    return(data)


def average_data(data, hybrid=False):
    '''Take replicate data files and average results across 
    each gene. Correct for differences in library size. Return 
    a dict of gene : data pairs.'''

    files = []
    averagedData = {}

    if len(data) == 0:
        return({})

    else:
        for dataFile in data:
            data_dict = data_to_dict(dataFile)
            files.append(data_dict)

        genes = list(data_dict.keys())
        genes.remove('librarySize')
        genes.sort()

        # For each gene, average across samples
        for gene in genes:
            refTotal = 0
            altTotal = 0
            refCounts = []
            altCounts = []
            pvals = []

            for i in files:
                librarySize = int(i['librarySize'])

                ref = int(i[gene][0])
                alt = int(i[gene][1])

                refCounts.append(ref)
                altCounts.append(alt)

                # Need to run binomial test on raw counts, so do it here
                if hybrid:
                    pvalue = binom_test(ref, ref + alt)
                    if pvalue < 0.05:
                        pvals.append(1)
                    else:
                        pvals.append(0)

                else:
                    pvals = [1]

                # Normalize for library size
                refTotal += (ref + 1) / librarySize
                altTotal += (alt + 1) / librarySize
            
            # Average expression values
            refAvg = refTotal / len(files)
            altAvg = altTotal / len(files)

            # Need these for minimum read cutoffs imposed later
            refMin = min(refCounts)
            altMin = min(altCounts)

            pvalAvg = sum(pvals) / len(pvals)

            averagedData[gene] = [refAvg, altAvg, refMin, altMin, pvalAvg]

        return(averagedData)


def merge_parents(parentA, parentB):
    '''Takes dict of expression data from each parental species 
    and returns one dict with genes as keys and log2(expression
    change) as values.'''

    merged = {}
    genesA = list(parentA.keys())
    genesB = list(parentB.keys())
    genes = list(set(genesA).intersection(genesB))
    genes.sort()

    for gene in genes:

        expA = parentA[gene][0]
        expB = parentB[gene][1]

        # Minimum read counts per sample
        if parentA[gene][2] > 5 and parentB[gene][3] > 5:
            value = log2(expA/expB)
            merged[gene] = value

        else:
            merged[gene] = 'NA'

    return(merged)


def filter_mapping_bias(genomeA, genomeB):
    '''Takes two dicts of hyrbid data mapped on each genome and
    returns a single dict of genomeA data with the biased genes
    removed. GenomeA is ideally the genome with the best assembly.'''

    unbiased = {}
    genesA = list(hybridA.keys())
    genesB = list(hybridB.keys())
    genes = list(set(genesA).intersection(genesB))
    genes.sort()

    for gene in genes:
        
        expAA = genomeA[gene][0]
        expAB = genomeA[gene][1]
        expBA = genomeB[gene][0]
        expBB = genomeB[gene][1]

        # Minimum cutoff of 20 reads per gene, half or more replicates significant
        if genomeA[gene][2] + genomeA[gene][3] > 20 and genomeA[gene][4] >= 0.5:

            valueA = log2(expAA/expAB)
            valueB = log2(expBA/expBB)

            # Cutoff for mapping bias (could play around with this)
            if abs(valueA-valueB) < 1.5:
                unbiased[gene] = valueA

            else: 
                unbiased[gene] = "NA"

        else: 
            unbiased[gene] = "NA"

    return(unbiased)


def filter_directionality(hybrid, parent):
    '''Filter genes for ones that match in directionality of effect 
    for parental and hybrid data.'''

    filtered = {}
    genes = list(hybrid.keys())
    genes.sort()

    for gene in genes:

        if hybrid[gene] == "NA":
            filtered[gene] = 'NA'

        elif parent[gene] == "NA":
            filtered[gene] = hybrid[gene]

        # Same directionality, keep it
        elif hybrid[gene]*parent[gene] > 0:
            filtered[gene] = hybrid[gene]
            
        else:
            filtered[gene] = "NA"

    return filtered


def write_dict(data, outFile):
    '''Take a dict of ASE data and write it to file.'''

    with open(outFile, 'w') as outFile:
        genes = list(data.keys())
        genes.sort()

        for gene in genes:
            if data[gene] != 'NA':
                outFile.write(gene + '\t' + str(data[gene]) + '\n')


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('-pa', '--parentA', 
        help="ASEr files from parentA expression data", 
        nargs="*")
    parser.add_argument('-pb', '--parentB', 
        help="ASEr files from parentB expression data", 
        nargs="*")
    parser.add_argument('-xa', '--hybridA', 
        help="ASEr files from hybrid expression data mapped on parentA genome", 
        nargs="*")
    parser.add_argument('-xb', '--hybridB', 
        help="ASEr files from hybrid expression data mapped on parentB genome", 
        nargs="*")
    parser.add_argument('-p', '--parentOut', 
        help="Output file for parental expression data")

    parser.add_argument('-x', '--hybridOut', 
        help="Output file for hybrid expression data")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()

    parentA = average_data(args.parentA)
    parentB = average_data(args.parentB)
    hybridA = average_data(args.hybridA, True)
    hybridB = average_data(args.hybridB, True)

    parentalExp = merge_parents(parentA, parentB)

    unbiasedHybrid = filter_mapping_bias(hybridA, hybridB)

    filteredHybrid = filter_directionality(unbiasedHybrid, parentalExp)

    write_dict(parentalExp, args.parentOut)
    write_dict(filteredHybrid, args.hybridOut)



