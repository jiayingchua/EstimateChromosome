## Take FASTA and return estimated number of chromosomes + their sizes

rm(list = ls())
graphics.off()

if(!require(seqinr)){
    install.packages("seqinr")
    library(seqinr)
}
args = commandArgs(trailingOnly = TRUE)

# argument 1 = fasta file path
# argument 2 = data frame output file  path
fastafile <- args[1]
outputtext <- args[2]

fastafile <- "C:\\Users\\JiaYing\\GP\\Anitra1161Contigs.fasta.gz"

# read fasta and return data frame with sequence id + lengths, assuming 'chromosomes' are atleast 1e+7 bp
fasta <- read.fasta(fastafile, seqtype = 'DNA')
lengths <- c()
header_length_table <- data.frame()
isChr <- c()
isChrlen <- c()
estChr <- data.frame()

lengths <- getLength(fasta)
header <- gsub("[^0-9A-Za-z]", "", getAnnot(fasta))
for (h in 1:length(header)) {
  header_length_table <-
    rbind(header_length_table, c(header[[h]], as.numeric(lengths[h])))
  if (lengths[h] > 1e+7) {
    isChr <- append(isChr, header[h])
    isChrlen <- append(isChrlen, lengths[h])
    estChr <-
      rbind(estChr, c(header[h], lengths[h]))
  }
}
colnames(estChr)[1] <- "seqid"
colnames(estChr) [2] <- "seqlen"

estChr <- estChr[order(estChr$seqid),]

write.table(estChr, file = outputtext, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
