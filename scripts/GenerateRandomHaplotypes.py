#!/usr/bin/python3

helpMessage = """
Generate "random" haplotypes for all bubble chains
of a Shasta phased assembly.

Each bubble chain generates two haplotypes obtained
by concatenating UR and PR contigs in Assembly-Phased.fasta
in the appropriate order.

Because UR contigs are not phased relative to each other,
this will generate switch errors.

Run this while in the assembly directory.
This uses as input PhasingRegions.csv and Assembly-Phased.fasta.
It generates output files Assembly-Random-Haplotype0.fasta
and Assembly-Random-Haplotype1.fasta

This script has no dependencies other than python3
and can invoked directly without any installation required.

"""

# Import what we need.
import argparse
import csv

# Make sure we have a --help option.
parser = argparse.ArgumentParser(description=helpMessage)
parser.parse_args()

# Read the bubble chains file.
csvFile = open('PhasingRegions.csv', 'r')
reader = csv.DictReader(csvFile)

bubbleChains = {}
for row in reader:
    bubbleChainId = int(row['Bubble chain id'])
    if not bubbleChainId in bubbleChains:
        bubbleChains[bubbleChainId] = []
    bubbleChains[bubbleChainId].append(row)
if not bubbleChains:
    print("No bubble chains were found."
        "Run this script from a Shasta phased assembly directory.")


# Read the Assembly-Phased.fasta file.
# Shasta writes each contig in a header line plus
# a single line containing sequence.
inputFastaFile = open('Assembly-Phased.fasta', 'r')
inputContigs = {}
while True:
    header = inputFastaFile.readline()
    if not header:
        break;
    if not header[0] == ">":
        raise RuntimeError("Invalid FASTA header: " + header)
    name = header[1:].split(" ")[0]
    sequence = inputFastaFile.readline().rstrip("\n")
    
    # We only want to keep it if the name begins with "UR," or "PR.". 
    if len(name) < 3:
        continue;
    prefix = name[0:3]
    if not (prefix == "UR." or prefix == "PR."):
        continue;
    
    inputContigs[name] = sequence
    

# Open the output files, one for each haplotype.
outputFileNames = [("Assembly-Random-Haplotype%i.fasta" % haplotypeId) for haplotypeId in range(2)]
outputFiles = [open(outputFileName, "w") for outputFileName in outputFileNames]


# Loop over bubble chains.
for bubbleChainId, bubbleChain in bubbleChains.items():
    print("Working on bubble chain %i of %i" % (bubbleChainId, len(bubbleChains)))
    
    # Check the bubble chain id.
    for x in bubbleChain:
        assert not x["Bubble chain id"] == bubbleChainId
        
    # Generate the two haplotypes for this bubble chain.
    for haplotypeId in range(2):
        sequence = ""
        for position in range(len(bubbleChain)):
            row = bubbleChain[position]
            if (row["Phased"]) == "No":
                name = "UR.%i.%i" % (bubbleChainId, position)
            else:
                component = int(row["Component"])                
                name = "PR.%i.%i.%i.%i" % (bubbleChainId, position, component, haplotypeId)
            assert name in inputContigs
            sequence += inputContigs[name]
        print("Bubble chain %i random haplotype %i has length %i" % (bubbleChainId, haplotypeId, len(sequence)))
        outputFiles[haplotypeId].write(">BC.%i.%i %i\n%s\n" % (bubbleChainId, haplotypeId, len(sequence), sequence))
            
print("Generation of random haplotypes is complete.")
print("These haplotypes can contain a switch error at each phased region.") 
print("Output is in %s and %s" % (outputFileNames[0], outputFileNames[1]))


