#! /usr/bin/python

import os
import sys
import gzip

args = sys.argv

PARSER = argparse.ArgumentParser(description='Program to extract exon/intron regions for target capture baitset designs.', usage'%(prog)s [options]', add_help=True)
PARSER.add_argument('geneList', metavar='genes', type=str, help='A comma delimited list of gene names from which to identify exon/intron regions. Note that these gene name should correspond to what exists in the annotation files.')
PARSER.add_argument('-o', '--output_dir', dest='outputDir', default='', help='Output directory to store output files. [default: %(default)s]')
PARSER.add_argument('-v', '--ensembl_ver', dest='ensemblVer', default='75', help='Ensembl annotation version to use. [default: %(default)s]')
PARSER.add_argument('-f', '--features', dest='features', default='exon,intron', help='Features to target in the gene (exon, intron). [default: %(default)s]')
PARSER.add_argument('-t', '--transcripts', dest='selectTrxs', default=None, help='List of transcript IDs to identify regions.')

geneNames = sys.argv[1].split(',')
regionType = sys.argv[2]
ensemblVer = '75'

class GTFFeature:
    def __init__(self, featureType, chrom, start, end, strand, meta):
        self.meta = meta
        self.featureType = featureType
        self.chrom = chrom.replace('chr', '')
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.geneId = meta['gene_id']
        self.exonId = None # meta['exon_id']
        self.exonNum = None #meta['exon_number']
        self.transcriptId = None #meta['transcript_id']
        self.geneName = meta['gene_name']
        self.transcriptName = None # meta['transcript_name']
        self.transcriptBioType = None # meta['transcript_biotype']
        self.source = meta['gene_source']
        self.geneBioType = meta['gene_biotype']
        self.set_values()

    def set_values(self):
        if self.featureType != "gene":
            self.transcriptId = self.meta['transcript_id']
            self.transcriptName = self.meta['transcript_name']
            self.transcriptBioType = self.meta['transcript_biotype']
            if self.featureType == 'exon':
                self.exonNum = self.meta['exon_number']
                self.exonId = self.meta['exon_id']
            elif self.featureType == 'CDS' or self.featureType == 'intron':
                self.exonNum = self.meta['exon_number']


class Gene:
    def __init__(self, geneId, gtfFeature):
        self.geneId = geneId
        self.features = {}

    def add_feature(self, gtfFeature):
        if gtfFeature.featureType not in self.features:
            self.features[gtfFeature.featureType] = []
        self.features[gtfFeature.featureType].append(gtfFeature)

    def get_regions(self, regionType):
        regions = None
        if regionType in self.features:
            regions = self.features[regionType]
        elif regionType == 'intron':
            regions = []
            for transcript in self.features['transcript']:
                transcriptRegions = []
                if transcript.transcriptBioType != 'protein_coding':
                    continue
                for exonRegion in self.get_transcript_regions(transcript.transcriptId):
                    transcriptRegions.append(exonRegion)
                tRegionsSorted = sorted(transcriptRegions, key=lambda x: x.start)
                for i in range(len(tRegionsSorted)-1):
                    intronStart = tRegionsSorted[i].end
                    intronEnd = tRegionsSorted[i+1].start
                    if tRegionsSorted[i].strand == '-':
                        tRegionsSorted[i].meta['exon_number'] = int(tRegionsSorted[i].meta['exon_number']) - 1
                    gf = GTFFeature("intron", tRegionsSorted[i].chrom, intronStart, intronEnd, tRegionsSorted[i].strand, tRegionsSorted[i].meta)
                    regions.append(gf)
        return regions

    def get_transcripts(self):
        return self.features['transcript']

    def get_transcript_regions(self, transcriptId):
        regions = self.get_regions('exon')
        transcriptRegions = []
        for region in regions:
            if region.transcriptId == transcriptId:
                transcriptRegions.append(region)
        return transcriptRegions


class GTF:
    def __init__(self, fn, ensemblVer):
        self.fn = fn
        self.ensVer = ensemblVer
        self.genes = {}
        self.parse()

    def parse(self):
        for line in open(self.fn):
            line = line.strip()
            if line.find('#!') > -1:
                continue
            self.process_line(line)

    def process_line(self, line):
        linesplit = line.split('\t')
        if linesplit[0].find('PATCH') > -1:
            return
        chrom, source, featureType, start, end, fill1, strand, fill2, meta = linesplit
        metaDict = self.parse_meta(meta)
        # print self.ensVer == '75'
        if self.ensVer == '75':
            metaDict['transcript_biotype'] = source
        # print line
        # print metaDict
        gf = GTFFeature(featureType, chrom, start, end, strand, metaDict)
        if featureType == "gene":
            if gf.geneId not in self.genes:
                self.genes[gf.geneName] = Gene(gf.geneId, gf)
        else:
            self.genes[gf.geneName].add_feature(gf)

    def parse_meta(self, metaValues):
        metaD = {}
        for value in metaValues.split('; '):
            # print value
            key, value = value.split(' ')
            value = value.rstrip('"').lstrip('"')
            metaD[key] = value
        return metaD

def get_feature_regions(gtf, gene, regionType):
    print gene
    regionOverlaps = []
    regions = gtf.genes[gene].get_regions(regionType)
    regionsSorted = sorted(regions, key=lambda x: x.start)
    for region in regionsSorted:
        # print region.chrom, region.start, region.end, region.geneName, region.transcriptId, region.transcriptBioType, region.geneBioType
        # Skip non-coding transcripts
    #            if region.transcriptBioType != 'protein_coding':
    #                continue
        overlap = False
        for regionOverlap in regionOverlaps:
            if region.chrom == regionOverlap[0]:
                if int(region.start) >= regionOverlap[1] and int(region.start) <= regionOverlap[2]:
                    # print region.chrom, region.start, region.end, 'overlaps with stored region', regionOverlap
                    overlap = True
                elif int(region.end) <= regionOverlap[2] and int(region.end) >= regionOverlap[1]:
                    # print region.chrom, region.start, region.end, 'overlaps with stored region', regionOverlap
                    overlap = True
                if overlap:
                    regionOverlap[3].append(region)
                    regionOverlap[1] = min(int(region.start), regionOverlap[1])
                    regionOverlap[2] = max(int(region.end), regionOverlap[2])
        if not overlap:
            # print 'No overlap, adding to list'
            regionOverlaps.append([region.chrom, int(region.start), int(region.end), [region]])

    regionIter = 1
    for regionOverlap in regionOverlaps:
#        print regionIter, regionOverlap[0], str(regionOverlap[1]) + "-" + str(regionOverlap[2])
        for region in regionOverlap[3]:
            print '\t'.join([str(x) for x in [regionIter, regionOverlap[0], str(regionOverlap[1]), str(regionOverlap[2]), region.chrom, region.start, region.end, region.transcriptId, region.exonNum]])
        regionIter += 1
    regionIter = 1
    for regionOverlap in regionOverlaps:
        tList = ','.join([x.transcriptId for x in regionOverlap[3]])
        exonNums = ','.join([str(x.exonNum) for x in regionOverlap[3]])
        print '\t'.join([str(x) for x in [regionIter, regionOverlap[0], str(regionOverlap[1]), str(regionOverlap[2]), tList, exonNums]])
        regionIter += 1

gtf_fn = '/data/ccgd/reference/human/gencode/GRCh37-p13/annotation/75/primary-assembly/Homo_sapiens.GRCh37.75.gtf.gz'
if ensemblVer == '82':
    gtf_fn = '/data/ccgd/reference/human/gencode/GRCh37-p13/annotation/82/primary-assembly/Homo_sapiens.GRCh37.82.chr.gtf.gz'
elif ensemblVer == '75':
    gtf_fn = '/data/ccgd/reference/human/gencode/GRCh37-p13/annotation/75/primary-assembly/Homo_sapiens.GRCh37.75.gtf.gz'

for gene in geneNames:
    gene_gtf = gene + "_tmp.gtf"
    cmd = 'zcat %s | grep -w %s > %s' % (gtf_fn, gene, gene_gtf)
    print cmd
    os.system(cmd)
    gtf = GTF(gene_gtf, ensemblVer)
    get_feature_regions(gtf, gene, regionType)
