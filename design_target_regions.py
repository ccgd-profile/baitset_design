#! /usr/bin/python

import os
import sys
import gzip
import logging
import argparse
import TargetDesigner.utils as utils
import TargetDesigner.transcript_annotation as transcript_annotation

args = sys.argv

PARSER = argparse.ArgumentParser(description='Program to extract exon/intron regions for target capture baitset designs.', usage='%(prog)s [options]', add_help=True)
PARSER.add_argument('geneList', metavar='genes', type=str, help='A comma delimited list of gene names from which to identify exon/intron regions. Note that these gene name should correspond to what exists in the annotation files.')
PARSER.add_argument('-o', '--output_dir', dest='outputDir', default='', help='Output directory to store output files. [default: %(default)s]')
PARSER.add_argument('-v', '--ensembl_ver', dest='ensemblVer', default='75', help='Ensembl annotation version to use. [default: %(default)s]')
PARSER.add_argument('-f', '--features', dest='features', default='exon,intron', help='Features to target in the gene (exon, intron). [default: %(default)s]')
PARSER.add_argument('-t', '--transcripts', dest='selectTrxs', default=None, help='List of transcript IDs to identify regions.')

pArgs = PARSER.parse_args()

print pArgs

genes = pArgs.geneList.split(',')
features = pArgs.features
ensemblVer = pArgs.ensemblVer

outDir = os.curdir()
if pArgs.outputDir is not None:
    outDir = pArgs.outputDir

if not os.path.isdir(outDir):
    os.mkdir(outDir)

featureTypes = pArgs.features.split(',')

gtf_fn = '/data/ccgd/reference/human/gencode/GRCh37-p13/annotation/75/primary-assembly/Homo_sapiens.GRCh37.75.gtf.gz'
if ensemblVer == '82':
    gtf_fn = '/data/ccgd/reference/human/gencode/GRCh37-p13/annotation/82/primary-assembly/Homo_sapiens.GRCh37.82.chr.gtf.gz'
elif ensemblVer == '75':
    gtf_fn = '/data/ccgd/reference/human/gencode/GRCh37-p13/annotation/75/primary-assembly/Homo_sapiens.GRCh37.75.gtf.gz'

for gene in genes:
    for featureType in featureTypes:
        if featureType == 'exon' or featureType == 'intron':
            continue
        geneDir = os.path.join(outDir, gene)
        geneOutFn = os.path.join(geneDir, gene + '_%s' % featureType)
        geneOutF = open(geneOutFn, 'w')
        if not os.path.isdir(geneDir):
            os.mkdir(geneDir)
        gene_gtf = os.path.join(geneDir, gene + "_tmp.gtf")
        cmd = 'zcat %s | grep -w %s > %s' % (gtf_fn, gene, gene_gtf)
        os.system(cmd)
        gtf = transcript_annotation.GTF(gene_gtf, ensemblVer)
        gtf.get_feature_regions(geneOutF, gene, regionType)
