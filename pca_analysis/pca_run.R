library("data.table")
library("seqinr")
library("mltools")
set.seed(20)
df2 <- fread("df1_aln.csv", stringsAsFactors=TRUE, header = T)
index <- seq(1,3665)

dfsamp <- as.data.frame(df2[index,])
countuniques <- as.numeric(lengths(lapply(dfsamp, unique)))
dfsamp <- dfsamp[,which(countuniques != 1)]
dfsamp <- as.data.table(dfsamp)
dfsamp$ID <- as.character(dfsamp$ID)
dfsamp <- one_hot(dfsamp)

fwrite(dfsamp, "dfsamp_aln_full.csv",quote = F)
pcinput <-dfsamp[,-1]
pcinput <- as.data.frame(pcinput)
pcoutput <- prcomp(pcinput,center = TRUE)
save(pcoutput, file = "pcoutput_full.RData")