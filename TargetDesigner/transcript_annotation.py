#! /usr/bin/python
# -*- coding: utf-8 -*-


class GTFFeature:
    """A class to contain all the elements of a particular transcript or transcript
    feature. This can be an exon, intron, UTR, or transcript.
    """

    def __init__(self, featureType, chrom, start, end, strand, meta):
        """Init function to set all the class variables. These values come directly
        from the GTF file.

        Args:
            featureType (string):
            chrom (string):
            start (int):
            end (int):
            strand (string):
            meta (string):

        Returns:
            None
        """

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
        """A function to parse the meta values from the GTF input values.

        It determines if the record passed in is a gene or not and sets the
        values accordingly.
        """

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
    """Gene class that contains different categories of features. These are encoded in the
    GTF file as transcript, exon, UTR.

    The features dictionary contains the featureType as the key and a list of the these
    GTFFeatures.

    When the exons for a gene are requested it returns the list in self.features['exon'].
    """

    def __init__(self, geneId, gtfFeature):
        """Initialize the Gene object with a gene ID and a features dictionary.

        Args:
            geneId (String):         The identifier for the gene as encoded in the GTF file meta field gene_id field.
            gtfFeature (GTFFeature): 
        """

        self.geneId = geneId
        self.features = {}

    def add_feature(self, gtfFeature):
        """
        """

        if gtfFeature.featureType not in self.features:
            self.features[gtfFeature.featureType] = []
        self.features[gtfFeature.featureType].append(gtfFeature)

    def get_regions(self, regionType):
        """
        """

        regions = None
        if regionType in self.features:
            regions = self.features[regionType]
        elif regionType == 'intron':
            regions = []
            for transcript in self.features['transcript']:
                transcriptRegions = []
                # Skip non-protein encoding transcripts.
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
        """Return all the transcripts for a gene. These will be coded as 
        GTFFeature objects.
        """

        return self.features['transcript']

    def get_transcript_regions(self, transcriptId):
        """Return exons (GTFFeature objects) for a specific transcript.
        """

        regions = self.get_regions('exon')
        transcriptRegions = []
        for region in regions:
            if region.transcriptId == transcriptId:
                transcriptRegions.append(region)
        return transcriptRegions


class GTF:
    """A class to organize the elements in the GTF file input.
    It organizes the annotation in the hierarchy that the GTF file contains.

    Gene:
        - features (exon, intron, transcript)
            - feature elements (ID, number, coordinates)
    """

    def __init__(self, fn, ensemblVer, trxList):
        """Instantiate the GTF object with the GTF filename, the Ensembl version
        and the transcript list.

        Parse the GTF file into gene objects and store all the gene features.
        """

        self.fn = fn
        self.ensVer = ensemblVer
        self.genes = {}
        self.trxList = trxList
        self.parse()

    def parse(self):
        """
        """

        for line in open(self.fn):
            line = line.strip()
            if line.find('#!') > -1:
                continue
            self.process_line(line)

    def process_line(self, line):
        """

        Args:
            line ():
        Returns:
            None
        """

        linesplit = line.split('\t')
        if linesplit[0].find('PATCH') > -1:
            return
        chrom, source, featureType, start, end, fill1, strand, fill2, meta = linesplit
        metaDict = self.parse_meta(meta)
        if self.ensVer == '75':
            metaDict['transcript_biotype'] = source
        gf = GTFFeature(featureType, chrom, start, end, strand, metaDict)
        if featureType == "gene":
            if gf.geneName not in self.genes:
                self.genes[gf.geneName] = Gene(gf.geneId, gf)
        else:
            addFeature = True
            if len(self.trxList) > 0:
                if gf.transcriptId not in self.trxList:
                    addFeature = False
            if addFeature:
                self.genes[gf.geneName].add_feature(gf)

    def parse_meta(self, metaValues):
        """
        """

        metaD = {}
        for value in metaValues.split('; '):
            # print value
            key, value = value.split(' ')
            value = value.rstrip('"').lstrip('"')
            metaD[key] = value
        return metaD

    def get_feature_regions(self, outF, gene, regionType, args):
        """

        Args:
            outF ():
            gene ():
            regionType ():
            args ():

        Returns:
        """

        regionOverlaps = []
        regions = self.genes[gene].get_regions(regionType)
        regionsSorted = sorted(regions, key=lambda x: x.start)
        for region in regionsSorted:
            overlap = False
            for regionOverlap in regionOverlaps:
                if region.chrom == regionOverlap['chrom']:
                    if int(region.start) >= regionOverlap['start'] and int(region.start) <= regionOverlap['end']:
                        overlap = True
                    elif int(region.end) <= regionOverlap['end'] and int(region.end) >= regionOverlap['start']:
                        overlap = True
                    if overlap:
                        regionOverlap['regionList'].append(region)
                        regionOverlap['start'] = min(int(region.start), regionOverlap['start'])
                        regionOverlap['end'] = max(int(region.end), regionOverlap['end'])
            if not overlap:
                regionOverlaps.append({'geneName': region.geneName,
                                       'chrom': region.chrom,
                                       'start': int(region.start),
                                       'end': int(region.end),
                                       'regionList': [region]}
                                      )

        regionIter = 1
        for regionOverlap in regionOverlaps:
            tList = ','.join([x.transcriptId for x in regionOverlap['regionList']])
            exonNums = ','.join([str(x.exonNum) for x in regionOverlap['regionList']])
            outF.write('\t'.join([str(x) for x in [regionOverlap['geneName'], regionIter, regionOverlap['chrom'], str(regionOverlap['start'] - args.upstreamBuffer), str(regionOverlap['end'] + args.downstreamBuffer), tList, exonNums]]) + '\n')
            regionIter += 1
