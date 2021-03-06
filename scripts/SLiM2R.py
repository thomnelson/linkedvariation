#!/usr/bin/env python
# 
# THIS SCRIPT TAKES SLiM OUTPUT FROM A GIVEN GENERATION AND
# CONVERTS IT TO A TABULAR FORMAT THAT IS EASILY PARSED IN R
#

import sys
import argparse

# FOR TESTING

parser = argparse.ArgumentParser(description='Convert SLiM output to tabular format for import to R')
parser.add_argument('-m','--mutations_file', required=True, help='File generated by a call to "outputMutations()".')
parser.add_argument('-g','--genome_files', required=True, nargs='+', help='List of files generated by calls to "outputSample()".')
parser.add_argument('-ng','--n_genomes', required=True, help='Number of genomes to analyze per subpopulation.')
args=parser.parse_args()

mutations_file = args.mutations_file
ngenomes       = int(args.n_genomes)
genome_files   = args.genome_files
nsubpops       = len(genome_files)

datadir  = "/Users/thom/Dropbox/1.Oregon/Data/active/simulations/forward/SLiM/simulations/migrationSelectionBalance/outputs/"
allneu   = "allneumu_gen100000.txt"
p1sample = "sample_p1_gen100000.txt"
p2sample = "sample_p2_gen100000.txt"

# READ MUTATIONS FILE TO GENERATE DICTIONARY OF MUTATIONS AND STATES

def readMutations(filepath, nsubpops, ngenomesPerSubpop):
    mudict   = {}
    mus      = open(filepath,'r')
    ngenomes = int(nsubpops) * int(ngenomesPerSubpop)
    for mu in mus:
        mu     = mu.strip().split(' ')
        gen          = mu[1]
        subpop       = mu[3]
        ID           = mu[4]
        if ID in mudict:
            continue
        type         = mu[5]
        position     = mu[6]
        s            = mu[7]
        h            = mu[8]
        originSubpop = mu[9]
        originGen    = mu[10]
        count        = mu[11]
        states       = []
        for i in range(ngenomes):
            states.append(0)
        mudict[ID] = {'type': type,'position': position,'s': s,'h': h, 
                          'originSubpop': originSubpop,'originGen': originGen,'states': states
                          }
    return(mudict)

x = readMutations(datadir+allneu, nsubpops, ngenomes)

### read genome samples, convert to dictionary
###  one element is key:value pairs for mutation ID -> mutation index
###  one is list of lists of mutation indices present in each genome

def readGenomes(filePath, ngenomes):
    genomes = open(filePath,'r')
    # get rid of header lines
    genomes.readline() ; genomes.readline()
    # read mutation info until get to 'Genomes' line of file
    mutations = {}
    mu_section = True
    while True:
        mutation = genomes.readline()
        if "Genomes" in mutation:
            break
        mutation = mutation.strip().split(' ')
        mutations[mutation[0]] = mutation[1]
    # genome section to the end
    genome = []
    for i in range(ngenomes):
        g = genomes.readline().strip().split(' ')
        genome.append(g)
    return({'mutations':mutations,'genomes':genome})

# p1genomes = readGenomes(filePath = datadir+p1sample, ngenomes = ngenomes)
p2genomes = readGenomes(filePath = datadir+p2sample, ngenomes = ngenomes)

def SLiM2tab(mutationsfile, ngenomes, list_of_sample_files):
    nsubpops = len(list_of_sample_files)
    output = readMutations(mutationsfile, nsubpops, ngenomes)
    # mutations object already has available slots for 
    # allelic states. Move through genome samples and
    # add derived states to mutations object in the 
    # appropriate locations
    for i in range(nsubpops):
        genomes = readGenomes(list_of_sample_files[i], ngenomes)
        mutDict = genomes['mutations']
        for j in range(ngenomes):
            genome = genomes['genomes'][j]
            # get correct list index for this genome
            outputIndex = j + (i * ngenomes)
            # convert mutation index to ID
            for k in range(2,len(genome)):
                # move through each mutation, ignoring leading genome identifiers
                mutIndex = genome[k]
                # get the mutation ID from its index and insert '1' in this 
                #  genome's slot to represent the derived state
                mutID    = mutDict[mutIndex]
                if mutID in output:
                    output[mutID]['states'][outputIndex] = 1
    outputTabular = []
    outputHeader  = "ID\ttype\tposition\ts\th\toriginSubpop\toriginGen"
    # complete header with genome identifiers
    for i in range(nsubpops):
        for j in range(ngenomes):
            genomeName = "\tp%s.genome%s"%(i+1,j+1)
            outputHeader = outputHeader+genomeName
    outputHeader = outputHeader+"\n"
    outputTabular.append(outputHeader)
    # loop through all mutations and convert info to a single string
    for mutationID in output:
        m = output[mutationID]
        type = m['type']
        position = m['position']
        s = m['s']
        h = m['h']
        originSubpop = m['originSubpop']
        originGen = m['originGen']
        lineout = "%s\t%s\t%s\t%s\t%s\t%s\t%s"%(mutationID,type,position,s,h,originSubpop,originGen)
        states = m['states']
        if sum(states) == 0:
            # none of the sampled chromosomes  have the derived allele. Skip it.
            continue
        for i in states:
            lineout = lineout+"\t"+str(i)
        lineout = lineout+"\n"
        outputTabular.append(lineout)
    return(outputTabular)

output = SLiM2tab(mutations_file, ngenomes, list_of_sample_files = genome_files)

for i in output:
    sys.stdout.write(i)
