## Take FASTA and return estimated number of chromosomes + their sizes

rm(list = ls())
graphics.off()

if(!require(seqinr)){
  install.packages("seqinr")
  library(seqinr)
}
if(!require(factoextra)){
  install.packages("factoextra")
  library(factoextra)
}
args = commandArgs(trailingOnly = TRUE)

# argument 1 = fasta file path
# argument 2 = data frame output file  path
fastafile <- args[1]
outputtext <- args[2]

# fastafile <- "C:\\Users\\JiaYing\\GP\\Anitra1161Contigs.fasta.gz"
# outputtext <- "C:\\Users\\JiaYing\\GP\\Anitra_OUT.csv"


# read fasta and return data frame with sequence id + lengths, assuming 'chromosomes' are atleast 1e+7 bp
fasta <- read.fasta(fastafile, seqtype = 'DNA')
lengths <- c()
header_length_table <- data.frame()
row_length_table <- data.frame()
isChr <- c()
isChrlen <- c()
estChr <- data.frame()

lengths <- getLength(fasta)
header <- gsub("[^0-9A-Za-z]", "", getAnnot(fasta))

for (h in 1:length(header)) {
  header_length_table <-
    rbind(header_length_table, c(header[[h]], as.numeric(lengths[h])))
  row_length_table <- 
    rbind(row_length_table, c(h, as.numeric(lengths[h])))
}

colnames(header_length_table)[1] <- "seqid"
colnames(header_length_table) [2] <- "seqlen"
colnames(row_length_table)[1] <- "seqid"
colnames(row_length_table) [2] <- "seqlen"
header_length_table$seqlen <- as.numeric(header_length_table$seqlen)
#plot(row_length_table)

if (nrow(header_length_table) < 20) {
  k <- 1
} else {
  k <- 3
}

set.seed(123)
km <- kmeans(row_length_table, k)
#fviz_cluster(km, data = row_length_table, palette = "jco")
cluster <- which.max(km[["centers"]][,2])

for (h in 1:length(header)) {
  if (km[["cluster"]][h] == cluster) {
    isChr <- append(isChr, header[h])
    isChrlen <- append(isChrlen, lengths[h])
    estChr <-
      rbind(estChr, c(header[h], lengths[h]))
  }
}

colnames(estChr)[1] <- "seqid"
colnames(estChr) [2] <- "seqlen"

estChr <- estChr[order(estChr$seqid),]
 
write.table(estChr, file = outputtext, sep = ",", col.names = FALSE, row.names = FALSE)
