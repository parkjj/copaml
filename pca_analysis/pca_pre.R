library("data.table")
library("ggplot2")
library("stringr")
library("RColorBrewer")
library("dplyr")
library("ggrepel")
library("seqinr")
library("mltools")
set.seed(20)

aln <- read.fasta("_out.200421030227825koE0shaIj6CwKgmM3Orpylsfnormal.pir")
metadata <- fread("metadata_reordered.txt")

alndf <-  data.frame(matrix(ncol = 0, nrow = 100840))#factor(vec, levels = c("a", "c", "g", "t", "n", "z"))
for(i in 1:length(aln)){
  print(i)
  vec <- as.character(aln[[i]])
  vec[which(vec == "-")] <- "z"
  vec[!which(vec %in% c("a", "c", "g", "t", "n", "z"))] <- "n"
  vec <- factor(vec, levels = c("a", "c", "g", "t", "n", "z"))
  alndf <- cbind(alndf, vec)
}
colnames(alndf) <- 1:length(aln)
df1 <- data.table( ID= metadata$`GenBank Accession`, t(alndf))
colnames(df1) <- c("ID", seq(1,100840))
fwrite(df1, "df1_aln.csv",quote = F)