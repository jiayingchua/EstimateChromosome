## Take FASTA and return estimated number of chromosomes + their sizes

rm(list = ls())
graphics.off()

require(seqinr)
#setwd("C:/Users/JiaYing/Group Project")

# # files to test
# fastafile <-
#   "C:/Users/JiaYing/Group Project/FASTA/Ridaeus_Ras1_v1.fasta.gz"
# fastafile <-
#   "C:/Users/JiaYing/Group Project/Rubus_argutus/Hillquist_genome_v1_purged_primary_contigs_HiC.fasta.gz"
# fastafile <-
#   "C:/Users/JiaYing/Group Project/FASTA/Rubus_occ_V3_10-12-17.fasta.gz"
# fastafile <-
#   "C:/Users/JiaYing/Group Project/FASTA/Rubus_chingii_Hu_fina.fa.gz"
# fastafile <-
#   "C:/Users/JiaYing/Group Project/FASTA/Burbank_genome_v1_purged_primary_contigs_HiC.fasta.gz"
# fastafile <-
#   "C:/Users/JiaYing/Group Project/FASTA/Rubus_Idaeus.genome.fa.gz"
# fastafile <-
#   "C:/Users/JiaYing/Group Project/FASTA/Rubus_occidentalis_masked_v1.1.fasta.gz"
# fastafile <-
#   "C:/Users/JiaYing/Group Project/FASTA/Rubus_occidentalis_v1.1.fasta.gz"
# fastafile <-
#   "C:/Users/JiaYing/Group Project/Ridaeus_Ras1_v1.0/Ridaeus_Ras1_scaffolds_yahs.fasta"

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

print(estChr)
