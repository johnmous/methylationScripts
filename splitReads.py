#!/usr/bin/python3.6
## Title: Group reads based on heterozygous SNPs locations. Each group of reads will be placed in a separate BAM file
## Author: I. Moustakas, i.moustakas@lumc.nl

import argparse
import pysam
import re
import os
import subprocess

## For all possible alleles in a genomic position, return a dictionary {allele => [Read Objects With this allele]}
def BaseToReads(samFile, chr, pos):
    pileups = samFile.pileup(chr, pos)

    baseToReadRecord = {}
    for pileupCol in pileups:
        for pileupRead in pileupCol.pileups:
            if not pileupRead.is_del and not pileupRead.is_refskip and pileupCol.pos == pos:
                aln = pileupRead.alignment
                base = aln.query_sequence[pileupRead.query_position]
                if base not in baseToReadRecord:
                    baseToReadRecord[base] = [aln]
                else:
                    baseToReadRecord[base].append(aln)
    return(baseToReadRecord)


if __name__ == '__main__':
    ### Parse input args
    parser = argparse.ArgumentParser(description='Group reads based on heterozygous SNPs locations. Each group of reads will be placed in a separate BAM file')
    parser.add_argument('--aln', type=str, help='Alignment file')
    parser.add_argument('--location', type=str, help='Allele location to split on, formatted as chr:pos')
    parser.add_argument('--thr', type=float, help='Threshold for allele frequency bellow which no BAM file will be created from this group of reads')
    parser.add_argument('--outpath', type=str, help='Path to place the output')
    args = parser.parse_args()

    ## Get command line arguments
    chr, pos = args.location.split(":")
    pos = int(pos) -1
    alignFile = args.aln
    outputPath = args.outpath
    thr = args.thr
    samFile = pysam.AlignmentFile(args.aln, 'rb')
    alleleToReadRecord = BaseToReads(samFile, chr, pos)
    if not os.path.exists(outputPath):
        os.makedirs(outputPath)

    allRecordCounts = 0
    recordsToWrite = {}
    for allele in alleleToReadRecord:
        allRecordCounts += len(alleleToReadRecord[allele])
    for allele in alleleToReadRecord:
        recordList = alleleToReadRecord[allele]
        if len(recordList) > thr*allRecordCounts:
            recordsToWrite[allele] = recordList

    ## Get the sample name out of the bam file name. Expects a Bismark output file
    str_search = re.search('.*/(.+)_bismark_bt2.sorted.bam', alignFile)
    sampleName = str_search.group(1)

    for allele in recordsToWrite:
        alleleSpecificFileID = "{0}_{1}_{2}_{3}".format(sampleName, chr, pos, allele)
        outputBam = "{0}/{1}_bismark_bt2.bam".format(outputPath, alleleSpecificFileID)
        outputSortedBam = "{0}/bamFiles/{1}_bismark_bt2.sorted.bam".format(outputPath, alleleSpecificFileID)
        if not os.path.exists(outputPath+"/bamFiles"):
            os.makedirs(outputPath+"/bamFiles")
        alleleSpecReadsBam = pysam.AlignmentFile(outputBam, "wb", template=samFile)
        for rcd in recordsToWrite[allele]:
            alleleSpecReadsBam.write(rcd)
        alleleSpecReadsBam.close()
        pysam.sort("-o", outputSortedBam, outputBam)
        pysam.index(outputSortedBam)
        subprocess.run(["/path/to/bismark-0.18.2/bismark_methylation_extractor",
                        "--output", outputPath, "--bedGraph", "--gzip", outputSortedBam])
        bismarkOutFile = "{0}/CpG_OB_{1}_bismark_bt2.sorted.txt.gz".format(outputPath, alleleSpecificFileID)
        subprocess.run(["/path/to/Rscript", "/path/to/methylationPattern.R",
                        bismarkOutFile, "/path/to/amplicon/table", outputPath+"/patterns" ])

    samFile.close()

## Example run 
## ./splitReads.py --aln /path/to/sample_bismark_bt2.sorted.bam --location chr19:56840683 --outpath /path/to/outputDir --thr 0.05




