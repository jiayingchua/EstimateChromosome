## This script takes a delta file as input and returns the estimated number
## of chromosomes of the query and the reference 

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
deltafile <- args[1]
outputtext <- args[2]

#deltafile <- "C:\\Users\\JiaYing\\GP\\Anitra1161Contigs.delta"
#outputtext <- "C:\\Users\\JiaYing\\GP\\Anitra1161Contigs_out.csv"


# reading the delta file and saving only the headers e.g. ">chr2 Rid2 33319654 33297897"
deltalines <- readLines(deltafile)

headers <- deltalines[grep("^>", deltalines)]
headers <- gsub(">", "", headers)
refid <- c()
reflen <- c()
for (i in 1:length(headers)) {
  refid <- append(refid,strsplit(headers[i], " ")[[1]][1])
  reflen <- append(reflen,strsplit(headers[i], " ")[[1]][3])
}

ref_table <- data.frame(refid, reflen)
ref_table[is.na(ref_table)] <- 0
new_ref_table <- c()
new_ref_table <- ref_table[!duplicated(ref_table),]


if (nrow(ref_table) < 20) {
  k <- 1
} else {
  k <- 3
}

row_length_table <- data.frame()
for (i in 1:nrow(new_ref_table)) {
  row_length_table <- rbind(row_length_table, c(i, as.numeric(new_ref_table$reflen[i])))
}

set.seed(123)
km <- kmeans(row_length_table, k)
#fviz_cluster(km, data = row_length_table, palette = "jco")
cluster <- which.max(km[["centers"]][,2])

isChr <- c()
isChrlen <- c()
estChr <- data.frame()
for (h in 1:nrow(row_length_table)) {
  if (km[["cluster"]][h] == cluster) {
    isChr <- append(isChr, new_ref_table$refid[h])
    isChrlen <- append(isChrlen, new_ref_table$reflen[h])
    estChr <-
      rbind(estChr, c(new_ref_table$refid[h], new_ref_table$reflen[h]))
  }
}

colnames(estChr)[1] <- "seqid"
colnames(estChr) [2] <- "seqlen"

estChr <- estChr[order(estChr$seqid),]

write.table(estChr, file = outputtext, sep = ",", col.names = FALSE, row.names = FALSE)