baitset_design
============

A general script to determine the region(s) of interest for designing target capture baits.

This script will take a list of gene names, and transcripts IDs, as input and identify the coordinates
of the exons or introns of the genes (and transcripts) of interest. It uses a GTF file with transcript
annotations for the coordinates of the exons to determine which regions to extract. When no transcript
ID(s) is specified for a particular gene, the script will output the union of overlapping exons in multiple
transcripts.


Installation
----------

- There is no installation process for this program. It is a python script that can be
downloaded and run on any system with [Python 2.7](https://www.python.org/download/releases/2.7)
installed.

Usage
---------

List the available command line parameters.
```
python <PATH_TO_DIRECTORY_SCRIPT_IS_LOCATED>/design_target_regions.py -h
```

Example to extract exon regions for genes in gene_list.txt using all ensembl 75 transcripts
for those genes.
```
python <PATH_TO_DIRECTORY_SCRIPT_IS_LOCATED>/design_target_regions.py -g gene_list.txt -o $PWD -f exon
```

Requirements
---------

### Gene input file
This file contains the genes of interest to design baits for the exons/intron regions. It can have 1 or 2 columns.
  1. Gene - HUGO gene symbol which contains the variant.
  2. Transcript ID - A comma-delimited list of transcript IDs corresponding to the gene. Only exon/intron regions will be
  extracted from these transcripts.

### GTF file
This is a annotation file containing all the genes and their associated transcripts, along with the corresponding coordinates 
for the genes, transcripts and their features. There are default files, but these assume the user is using the Eris compute
cluster (in which the default files are located). There is an option (-a, --annotation) to input a different GTF file. Note that
the annotation file is assumed to be gzipped (.gz).

Parameters
-------------
| Parameter | Description | Default |
| ---------- | ----------- | ------- |
| -h, --help | Show the help message and exit | NA |
| -g, --genes | A file containing HUGO gene symbols of interest. The gene name need to match the gene names in the GTF annotation file. | NA (Required parameter) |
| -o, --output_dir | Output directory to store output files. | '.' |
| -v, --ensembl_ver | Ensembl annotation version to use. This assumes user is on the Eris computer cluster and the file paths are hardcoded. Use -a option for user-specified GTF file. | 75 |
| -f, --features | Features to target in the gene (exon, intron). | exon,intron |
| -u, --upstream_buffer | Number of base pairs that should be added upstream the regions of interest. | 0 |
| -d, --downstream_buffer | Number of base pairs that should be added downstream the regions of interest. | 0 |
| -a, --annotation | A gzipped GTF file containing transcript annotations for determining the regions to extract. | None |

Output files and formats
-----------

### Output file
When the program completes, there should be output files in gene-specific directories within the specified output_dir. There will be files for each featureType specified (i.e., exon and intron). The 
7 columns contain:
  1. gene name - The HUGO gene symbol that was input.
  2. region index - The index of the region to be baited.
  3. chrom - The chromosome of the region.
  4. start - The region start position.
  5. end - The region end position.
  6. transcript IDs - A comma-delimited list of transcript IDs containing regions covered by the region.
  7. feature index - A comma-delimited list of the feature indices. These are relative to the other features in the transcript. If extracting exons, this number would be the exon number.
