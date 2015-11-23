#! /usr/bin/python

import os
import sys
import gzip
import logging
import argparse
import TargetDesigner.utils as utils
import TargetDesigner.transcript_annotation as transcript_annotation

# ------------Setup parser -----------------------
PARSER = argparse.ArgumentParser(description='Program to extract exon/intron regions for target capture baitset designs.',
                                 usage='%(prog)s [options]', add_help=True)
PARSER.add_argument('-g', '--genes', help='A file containing a list of gene names from which to identify exon/intron regions. Note that \
                     these gene name should correspond to what exists in the annotation files. A second column can be input with specific \
                     transcript IDs for each gene.', required=True)
PARSER.add_argument('-o', '--output_dir', dest='outputDir', default='', help='Output directory to store output files. [default: %(default)s]')
PARSER.add_argument('-v', '--ensembl_ver', dest='ensemblVer', default='75', help='Ensembl annotation version to use. [default: %(default)s]')
PARSER.add_argument('-f', '--features', dest='features', default='exon,intron', help='Features to target in the gene (exon, intron). [default: %(default)s]')
PARSER.add_argument('-t', '--transcripts', dest='selectTrxs', default=None, help='List of transcript IDs to identify regions. [default: %(default)s]')
PARSER.add_argument('-u', '--upstream_buffer', dest='upstreamBuffer', default=0, type=int, help='Number of base pairs that should be added upstream the regions of interest. [default: %(default)s]')
PARSER.add_argument('-d', '--downstream_buffer', dest='downstreamBuffer', default=0, type=int, help='Number of base pairs that should be added downstream the regions of interest. [default: %(default)s]')
PARSER.add_argument('-a', '--annotation', dest='gtfFn', default=None, help='A gzipped GTF file containing transcript annotations for determining the regions to extract.')

pArgs = PARSER.parse_args()

# Parse the input file with gene names and possibly associated
# transcript IDs for each gene. Store these in a dictionary with the
# gene name as the key and the transcript list as the value.
inputGenes = utils.parse_input_gene_file(pArgs.genes)

features = pArgs.features
ensemblVer = pArgs.ensemblVer

# Create the specified output directory for all the files.
# If no output directory is specified it will use the current directory.
outDir = os.curdir
if pArgs.outputDir is not None:
    outDir = pArgs.outputDir

if not os.path.isdir(outDir):
    os.mkdir(outDir)

# Allow the specification of exons and introns (i.e., exon,intron).
featureTypes = pArgs.features.split(',')
# Check to make sure the featureType is exon or intron, otherwise
# exit with error message.
for featureType in featureTypes:
    if featureType not in ['exon', 'intron']:
        print 'An unsupported feature type was used, %s, exiting.' % featureType
        sys.exit(1)

# Use gzip GTF formatted files containing transcript annotations.
# These are hard-coded for use on the Eris cluster. If another system is
# used, then use the -a, --annotation option to pass in another file.
gtfFn = '/data/ccgd/reference/human/gencode/GRCh37-p13/annotation/75/primary-assembly/Homo_sapiens.GRCh37.75.gtf.gz'
if pArgs.annotation is None:
    if ensemblVer == '82':
        gtfFn = '/data/ccgd/reference/human/gencode/GRCh37-p13/annotation/82/primary-assembly/Homo_sapiens.GRCh37.82.chr.gtf.gz'
    elif ensemblVer == '75':
        gtfFn = '/data/ccgd/reference/human/gencode/GRCh37-p13/annotation/75/primary-assembly/Homo_sapiens.GRCh37.75.gtf.gz'
else:
    gtfFn = pArgs.gtfFn

print len(inputGenes.keys()), "genes in input list"

for geneName in inputGenes.keys():
    trxList = inputGenes[geneName]
    for featureType in featureTypes:
        geneDir = os.path.join(outDir, geneName)
        if not os.path.isdir(geneDir):
            os.mkdir(geneDir)
        geneOutFn = os.path.join(geneDir, geneName + '_%s' % featureType + '.txt')
        geneOutFile = open(geneOutFn, 'w')
        geneGTFFn = os.path.join(geneDir, geneName + "_tmp.gtf")
        cmd = 'zcat %s | grep -w %s > %s' % (gtfFn, geneName, geneGTFFn)
        os.system(cmd)
        gtf = transcript_annotation.GTF(geneGTFFn, ensemblVer, trxList)
        gtf.get_feature_regions(geneOutFile, geneName, featureType, pArgs)
