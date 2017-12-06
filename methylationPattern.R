## Author: I. Moustakas  
## Title: Get the methylation patterns and count them
## Usage: methylationPattern.R bismarkCpG outputDir sampleName
library("reshape2")
library("stringr")

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3 ) {
  stop("CpG file from bismark (CpG_OB_*), amplicons table and output direcotry must be supplied", call.=FALSE)
} 

bismarkCpG_OB = args[1]
ampliconsFile = args[2]
outDir = args[3]

dir.create(outDir, showWarnings = FALSE)
regionReportFile = paste0(outDir, "/regionReport.txt")
file.create(regionReportFile, showWarnings = FALSE)

## Get sample name out of the bismakr file name 
sampleName = str_match(bismarkCpG_OB, ".+/CpG_OB_(.+)_bismark.+")[1,2]
line <- sprintf("=========\nThis concerns sample: %s\n=========",sampleName )
regionReportFile = paste0(outDir, "/", sampleName, "_regionReport.txt")
file.create(regionReportFile, showWarnings = FALSE)
write(line, file=regionReportFile, append=TRUE)

## Load data. Forward and reverse strand are in two separate files (OB and OT). Combine them in one df
## If file does not exist, create empty DF 
emptyDF <- data.frame(V1 = factor(),
                      V2 = factor(),
                      V3 = factor(),
                      V4 = integer(),
                      V5 = factor())
if (file.exists(bismarkCpG_OB)) {
  methylationOB <- read.table(bismarkCpG_OB, skip = 1, sep="\t")
} else methylationOB <- emptyDF
bismarkCpG_OT <- gsub("CpG_OB_", "CpG_OT_", bismarkCpG_OB)
if (file.exists(bismarkCpG_OT)) {
  methylationOT <- read.table(bismarkCpG_OT, skip = 1, sep="\t")
}else methylationOT <- emptyDF
methylation <- rbind(methylationOB, methylationOT)
colnames(methylation) <- c("Read", "MethylStatus", "Chr", "Pos", "Zz")

## Load amplicon positions
amplicons <- read.table(ampliconsFile, header = TRUE)

SumPos <- function(methCounts, pattern){
  sapply(methCounts, function(c){
    if (is.na(c[pattern])) 0
    else c[pattern]
  })
}

SumCountPos <- function(sumMethPos, countPatterns) {
  sapply(sumMethPos, function(sum){
  sum(countPatterns[sumMethPos==sum])
}) 
}

## Go through the amplicons and extract from the methylation table
result <- apply(amplicons, 1, function(amplicon) {
  name = amplicon["Name"]
  chr = amplicon["Chr"]
  start = as.integer(amplicon["start"])
  end = as.integer(amplicon["end"])
  strand = amplicon["strand"]

  
  ampliconMethyl <- methylation[methylation$Chr == chr,] 
  ampliconMethyl <- ampliconMethyl[ampliconMethyl$Pos>=start & ampliconMethyl$Pos<=end, ]
  
  ## If there are records, proceed
  if (nrow(ampliconMethyl) > 0) {
    ampliconMethyl <- ampliconMethyl[ ,c(1,4,2)]
    methylPos <- ampliconMethyl$Pos
    countPos <- table(methylPos)
    sumRecords <- length(unique(ampliconMethyl[,1]))
    highCountPos <- names(countPos[countPos > 0.11*sumRecords])
    ampliconMethyl <- ampliconMethyl[ampliconMethyl$Pos %in% highCountPos, ]
    methylPattern <- dcast(ampliconMethyl, Read ~ Pos )
    positions <- as.character(paste0(sort(unique(ampliconMethyl$Pos)), ", "))
    numberPositions <- length(positions)
    line <- sprintf("\nIn Amplicon: %s on %s:%d-%d there were %d CpG positions detected. These are:" , name, chr, start, end, numberPositions) 
    write(line, file=regionReportFile, append=TRUE)
    write(positions, file=regionReportFile, append=TRUE)
    methylPattern[is.na(methylPattern)] <- "*"
    methylPattern <- methylPattern[-1]
    countPatterns <- table(do.call(paste0, methylPattern))
    listOfPatterns <- strsplit(names(countPatterns), NULL)
    patterns <- as.data.frame(do.call(rbind, listOfPatterns))
    
    methCounts<- apply(patterns,1, function(pattern){
      table(pattern)
    })

    colnames(patterns) <- colnames(methylPattern)
    patterns$counts <- as.vector(countPatterns)
    
    sumMethPos<- SumPos(methCounts, "+")
    patterns$sumMethPos <- sumMethPos
    sumCountsMethPos <- SumCountPos(sumMethPos,countPatterns)
    patterns$sumCountsMethPos <- sumCountsMethPos
    patterns$pcntSumCountsMethPos <- round(sumCountsMethPos/sumRecords, 3)
    
    sumUnMethPos <- SumPos(methCounts, "-")
    patterns$sumUnMethPos <- sumUnMethPos
    sumCountsUnMethPos <- SumCountPos(sumUnMethPos,countPatterns)
    patterns$sumCountsUnMethPos <- sumCountsUnMethPos
    patterns$pcntSumCountsUnMethPos <- round(sumCountsUnMethPos/sumRecords, 3)
    
    sumUnknownPos <- SumPos(methCounts, "*")
    patterns$sumUnknownPos <- sumUnknownPos
    sumCountsUnknownPos <- SumCountPos(sumUnknownPos,countPatterns)
    patterns$sumCountsUnknownPos <- sumCountsUnknownPos
    patterns$pcntSumCountsUnknownPos <- round(sumCountsUnknownPos/sumRecords, 3)
    patternFile <- sprintf("%s/%s_%s_methylation.tsv", outDir, sampleName, name)
    write.table(patterns, file=patternFile, quote = F, sep = "\t", row.names = F )
    
  } else {
    line <- sprintf("\nAmplicon %s Not Found", name)
    write(line, file=regionReportFile, append=TRUE)
    
  } 
})

