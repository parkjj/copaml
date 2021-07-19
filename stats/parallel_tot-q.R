library(stringr)
library(slider)
library(foreach)
library(parallel)
library(doParallel)

cl <- parallel::makeCluster(6)
doParallel::registerDoParallel(cl)

pooledscores_human <- read.csv("pooledscores_human.csv")
pooledscores_intersect <- read.csv("pooledscores_intersect.csv")
pooledscores_patho <- read.csv("pooledscores_patho.csv")

pooledscores <- cbind(pooledscores_human, pooledscores_intersect, pooledscores_patho)

cumentropy <- read.csv("../classprint_cumentropy_for_6bp_window.csv", header=F)
cumentropy_sub <- cumentropy[seq(1,dim(pooledscores)[1]),]
cumentropy_sub <- as.numeric(cumentropy_sub)
null_index <- which(cumentropy_sub <= min(cumentropy_sub))

null_distr <- as.numeric(unlist(pooledscores[null_index,]))
slideindex <- slide(seq(1, 100840), ~.x, .before = 5)

pvalues <- NULL

pvalues <- foreach(i = seq(1, 100840), .combine = 'c') %dopar% {
  play <- as.numeric(unlist(pooledscores[unlist(slideindex[i]),]))
  if(any(is.na(play))){
    play <- play[- which(is.na(play))]
  }
  p <- wilcox.test(play, null_distr, alternative = "two.sided")
  p$p.value
}
pvalues <- p.adjust(pvalues, method = "BH")

write.table(pvalues, "tot-q.csv", quote = F, row.names = F, col.names = F, sep = ",")
