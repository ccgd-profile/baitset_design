#! /usr/bin/python


def parse_input_gene_file(inputGeneFile):
    """
    """

    geneDict = {}
    for line in open(inputGeneFile, 'r'):
        line = line.strip()
        linesplit = line.split('\t')
        trxList = []
        if len(linesplit) > 1:
            trxList = linesplit[1].split(',')
        geneDict[linesplit[0]] = trxList

    return geneDict
