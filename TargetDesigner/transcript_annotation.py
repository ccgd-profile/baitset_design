#! /usr/bin/python
# -*- coding: utf-8 -*-


class GTFFeature:
    def __init__(self, featureType, chrom, start, end, strand, meta):
        self.meta = meta
        self.featureType = featureType
        self.chrom = chrom.replace('chr', '')
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.geneId = meta['gene_id']
        self.exonId = None
        self.exonNum = None
        self.transcriptId = None
        self.geneName = meta['gene_name']
        self.transcriptName = None
        self.transcriptBioType = None
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
        print 'Adding feature', gtfFeature.geneId, gtfFeature.transcriptId
        if gtfFeature.featureType not in self.features:
            self.features[gtfFeature.featureType] = []
        self.features[gtfFeature.featureType].append(gtfFeature)

    def get_regions(self, regionType):
        regions = None
        print 'Gene.get_regions()', self.features
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
                for i in range(len(tRegionsSorted) - 1):
                    intronStart = tRegionsSorted[i].end
                    intronEnd = tRegionsSorted[i + 1].start
                    if tRegionsSorted[i].strand == '-':
                        tRegionsSorted[i].meta['exon_number'] = int(tRegionsSorted[i].meta['exon_number']) - 1
                    gf = GTFFeature("intron", tRegionsSorted[i].chrom, intronStart, intronEnd, tRegionsSorted[i].strand, tRegionsSorted[i].meta)
                    regions.append(gf)
        if regions is None:
            print 'No regions found for', self.geneId
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
    def __init__(self, fn, ensemblVer, trxList):
        self.fn = fn
        self.ensVer = ensemblVer
        self.genes = {}
        self.trxList = trxList
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
            print 'Adding gene', line
            if gf.geneName not in self.genes:
                self.genes[gf.geneName] = Gene(gf.geneId, gf)
        else:
            addFeature = True
            if len(self.trxList) > 0:
                if gf.transcriptId not in self.trxList:
                    addFeature = False
            if addFeature:
                print 'Adding adding transcript element', line
                print 'Gene name', gf.geneName
                self.genes[gf.geneName].add_feature(gf)

    def parse_meta(self, metaValues):
        metaD = {}
        for value in metaValues.split('; '):
            # print value
            key, value = value.split(' ')
            value = value.rstrip('"').lstrip('"')
            metaD[key] = value
        return metaD

    def get_feature_regions(self, outF, gene, regionType):
        regionOverlaps = []
        regions = self.genes[gene].get_regions(regionType)
        regionsSorted = sorted(regions, key=lambda x: x.start)
        for region in regionsSorted:
            # print region.chrom, region.start, region.end, region.geneName, region.transcriptId, region.transcriptBioType, region.geneBioType
            # Skip non-coding transcripts
                    # if region.transcriptBioType != 'protein_coding':
                        # continue
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
            # print regionIter, regionOverlap[0], str(regionOverlap[1]) + "-" + str(regionOverlap[2])
            for region in regionOverlap[3]:
                print '\t'.join([str(x) for x in [region.geneName, regionIter, regionOverlap[0], str(regionOverlap[1]), str(regionOverlap[2]), region.chrom, region.start, region.end, region.transcriptId, region.exonNum]])
            regionIter += 1
        regionIter = 1
        for regionOverlap in regionOverlaps:
            tList = ','.join([x.transcriptId for x in regionOverlap[3]])
            exonNums = ','.join([str(x.exonNum) for x in regionOverlap[3]])
            outF.write('\t'.join([str(x) for x in [regionIter, regionOverlap[0], str(regionOverlap[1]), str(regionOverlap[2]), tList, exonNums]]))
            regionIter += 1
