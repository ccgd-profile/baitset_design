#! /usr/bin/python

import os
import sys
import gzip
import logging
import argparse
import TargetDesigner.utils as utils

args = sys.argv

PARSER = argparse.ArgumentParser(description='Program to extract exon/intron regions for target capture baitset designs.', usage='%(prog)s [options]', add_help=True)
PARSER.add_argument('geneList', metavar='genes', type=str, help='A comma delimited list of gene names from which to identify exon/intron regions. Note that these gene name should correspond to what exists in the annotation files.')
PARSER.add_argument('-o', '--output_dir', dest='outputDir', default='', help='Output directory to store output files. [default: %(default)s]')
PARSER.add_argument('-v', '--ensembl_ver', dest='ensemblVer', default='75', help='Ensembl annotation version to use. [default: %(default)s]')
PARSER.add_argument('-f', '--features', dest='features', default='exon,intron', help='Features to target in the gene (exon, intron). [default: %(default)s]')
PARSER.add_argument('-t', '--transcripts', dest='selectTrxs', default=None, help='List of transcript IDs to identify regions.')

pArgs = PARSER.parse_args()

print pArgs

genes = pArgs.geneList
features = pArgs.features
ensemblVer = pArgs.ensemblVer

sys.exit()

# gtf_fn = '/data/ccgd/reference/human/gencode/GRCh37-p13/annotation/75/primary-assembly/Homo_sapiens.GRCh37.75.gtf.gz'
# if ensemblVer == '82':
#     gtf_fn = '/data/ccgd/reference/human/gencode/GRCh37-p13/annotation/82/primary-assembly/Homo_sapiens.GRCh37.82.chr.gtf.gz'
# elif ensemblVer == '75':
#     gtf_fn = '/data/ccgd/reference/human/gencode/GRCh37-p13/annotation/75/primary-assembly/Homo_sapiens.GRCh37.75.gtf.gz'

# for gene in geneNames:
#     gene_gtf = gene + "_tmp.gtf"
#     cmd = 'zcat %s | grep -w %s > %s' % (gtf_fn, gene, gene_gtf)
#     print cmd
#     os.system(cmd)
#     gtf = GTF(gene_gtf, ensemblVer)
#     get_feature_regions(gtf, gene, regionType)


# def get_feature_regions(gtf, gene, regionType):
#     print gene
#     regionOverlaps = []
#     regions = gtf.genes[gene].get_regions(regionType)
#     regionsSorted = sorted(regions, key=lambda x: x.start)
#     for region in regionsSorted:
#         # print region.chrom, region.start, region.end, region.geneName, region.transcriptId, region.transcriptBioType, region.geneBioType
#         # Skip non-coding transcripts
#     #            if region.transcriptBioType != 'protein_coding':
#     #                continue
#         overlap = False
#         for regionOverlap in regionOverlaps:
#             if region.chrom == regionOverlap[0]:
#                 if int(region.start) >= regionOverlap[1] and int(region.start) <= regionOverlap[2]:
#                     # print region.chrom, region.start, region.end, 'overlaps with stored region', regionOverlap
#                     overlap = True
#                 elif int(region.end) <= regionOverlap[2] and int(region.end) >= regionOverlap[1]:
#                     # print region.chrom, region.start, region.end, 'overlaps with stored region', regionOverlap
#                     overlap = True
#                 if overlap:
#                     regionOverlap[3].append(region)
#                     regionOverlap[1] = min(int(region.start), regionOverlap[1])
#                     regionOverlap[2] = max(int(region.end), regionOverlap[2])
#         if not overlap:
#             # print 'No overlap, adding to list'
#             regionOverlaps.append([region.chrom, int(region.start), int(region.end), [region]])

#     regionIter = 1
#     for regionOverlap in regionOverlaps:
# #        print regionIter, regionOverlap[0], str(regionOverlap[1]) + "-" + str(regionOverlap[2])
#         for region in regionOverlap[3]:
#             print '\t'.join([str(x) for x in [regionIter, regionOverlap[0], str(regionOverlap[1]), str(regionOverlap[2]), region.chrom, region.start, region.end, region.transcriptId, region.exonNum]])
#         regionIter += 1
#     regionIter = 1
#     for regionOverlap in regionOverlaps:
#         tList = ','.join([x.transcriptId for x in regionOverlap[3]])
#         exonNums = ','.join([str(x.exonNum) for x in regionOverlap[3]])
#         print '\t'.join([str(x) for x in [regionIter, regionOverlap[0], str(regionOverlap[1]), str(regionOverlap[2]), tList, exonNums]])
#         regionIter += 1