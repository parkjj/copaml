library("data.table")
library("ggplot2")
library("stringr")
library("RColorBrewer")
library("dplyr")
library("ggrepel")
library("seqinr")
library("mltools")
set.seed(20)

metadata <- fread("metadata_reordered.txt")
df4 <- fread("dfsamp_aln_full.csv")
load("pcoutput_full.RData")

summary(pcoutput)
play <- summary(pcoutput)$importance
pc <- pcoutput$x
pc <- data.frame(df4$ID, pc[,1], pc[,2],  pc[,3], pc[,4])
colnames(pc) <- c("sample", "PC1", "PC2", "PC3", "PC4")

# host
pc$label <- metadata$Host[match(pc$sample, metadata$`GenBank Accession`)]
hostnames <- read.table("host.txt")
hostnames <- hostnames$V1
pc$label[- which(pc$label %in% hostnames)] <- "zz_other"

qcpc <- ggplot(data=pc, aes(x=PC1, y=PC2)) + geom_point( aes(color=label),  size = 5,pch=1, stroke = 0.4) # + scale_fill_manual(values=c("#9399ff", "#a9fffd", "#46cdcf",  "#f67280", "#c06c84", "#6c5b7b", "#355c7d"))# + #+ ggtitle("Boxplot of the sgRNA representation")
qcpc <- qcpc + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size = 0.8), axis.text=element_text(size=15), axis.title=element_text(size=20), axis.line = element_line(size = 0),  legend.title=element_text(size=0), legend.text=element_text(size=20), plot.title = element_text(size=22))
qcpc <- qcpc + labs(x = "PC1 (24%)", y = "PC2 (19%)") 
qcpc <- qcpc + scale_shape_manual(values=c(21, 2, 5)) #+ ggtitle("PCA Su CD4 dataset, 3-mers")
pdf("pca_full_host.pdf", height = 5, width = 7.5)
qcpc
dev.off()

qcpc <- ggplot(data=pc, aes(x=PC3, y=PC4)) + geom_point( aes(color=label),  size = 5,pch=1, stroke = 0.4) # + scale_fill_manual(values=c("#9399ff", "#a9fffd", "#46cdcf",  "#f67280", "#c06c84", "#6c5b7b", "#355c7d"))# + #+ ggtitle("Boxplot of the sgRNA representation")
qcpc <- qcpc + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size = 0.8), axis.text=element_text(size=15), axis.title=element_text(size=20), axis.line = element_line(size = 0),  legend.title=element_text(size=0), legend.text=element_text(size=20), plot.title = element_text(size=22))
qcpc <- qcpc + labs(x = "PC3 (13.7%)", y = "PC4 (10.6%)") 
qcpc <- qcpc + scale_shape_manual(values=c(21, 2, 5)) #+ ggtitle("PCA Su CD4 dataset, 3-mers")
pdf("pca_full_host_set2.pdf", height = 5, width = 7.5)
qcpc
dev.off()


# genus
pc$label <- metadata$`Virus Genus`[match(pc$sample, metadata$`GenBank Accession`)]
hostnames <- read.table("genus.txt")
hostnames <- hostnames$V1
pc$label[- which(pc$label %in% hostnames)] <- "zz_other"

qcpc <- ggplot(data=pc, aes(x=PC1, y=PC2)) + geom_point( aes(color=label),  size = 5,pch=1, stroke = 0.4) # + scale_fill_manual(values=c("#9399ff", "#a9fffd", "#46cdcf",  "#f67280", "#c06c84", "#6c5b7b", "#355c7d"))# + #+ ggtitle("Boxplot of the sgRNA representation")
qcpc <- qcpc + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size = 0.8), axis.text=element_text(size=15), axis.title=element_text(size=20), axis.line = element_line(size = 0),  legend.title=element_text(size=0), legend.text=element_text(size=20), plot.title = element_text(size=22))
qcpc <- qcpc + labs(x = "PC1 (24%)", y = "PC2 (19%)") 
qcpc <- qcpc + scale_shape_manual(values=c(21, 2, 5)) #+ ggtitle("PCA Su CD4 dataset, 3-mers")
pdf("pca_full_genus.pdf", height = 5, width = 9)
qcpc
dev.off()

qcpc <- ggplot(data=pc, aes(x=PC3, y=PC4)) + geom_point( aes(color=label),  size = 5,pch=1, stroke = 0.4) # + scale_fill_manual(values=c("#9399ff", "#a9fffd", "#46cdcf",  "#f67280", "#c06c84", "#6c5b7b", "#355c7d"))# + #+ ggtitle("Boxplot of the sgRNA representation")
qcpc <- qcpc + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size = 0.8), axis.text=element_text(size=15), axis.title=element_text(size=20), axis.line = element_line(size = 0),  legend.title=element_text(size=0), legend.text=element_text(size=20), plot.title = element_text(size=22))
qcpc <- qcpc + labs(x = "PC3 (13.7%)", y = "PC4 (10.6%)") 
qcpc <- qcpc + scale_shape_manual(values=c(21, 2, 5)) #+ ggtitle("PCA Su CD4 dataset, 3-mers")
pdf("pca_full_genus_set2.pdf", height = 5, width = 9)
qcpc
dev.off()

# species
pc$label <- metadata$`Virus Species`[match(pc$sample, metadata$`GenBank Accession`)]
hostnames <- read.csv("species.txt", header = F)
hostnames <- hostnames$V1
pc$label[- which(pc$label %in% hostnames)] <- "zz_other"

qcpc <- ggplot(data=pc, aes(x=PC1, y=PC2)) + geom_point( aes(color=label),  size = 5,pch=1, stroke = 0.4) # + scale_fill_manual(values=c("#9399ff", "#a9fffd", "#46cdcf",  "#f67280", "#c06c84", "#6c5b7b", "#355c7d"))# + #+ ggtitle("Boxplot of the sgRNA representation")
qcpc <- qcpc + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size = 0.8), axis.text=element_text(size=15), axis.title=element_text(size=20), axis.line = element_line(size = 0),  legend.title=element_text(size=0), legend.text=element_text(size=20), plot.title = element_text(size=22))
qcpc <- qcpc + labs(x = "PC1 (24%)", y = "PC2 (19%)") 
qcpc <- qcpc + scale_shape_manual(values=c(21, 2, 5)) #+ ggtitle("PCA Su CD4 dataset, 3-mers")
pdf("pca_full_species.pdf", height = 5, width = 13)
qcpc
dev.off()

qcpc <- ggplot(data=pc, aes(x=PC3, y=PC4)) + geom_point( aes(color=label),  size = 5,pch=1, stroke = 0.4) # + scale_fill_manual(values=c("#9399ff", "#a9fffd", "#46cdcf",  "#f67280", "#c06c84", "#6c5b7b", "#355c7d"))# + #+ ggtitle("Boxplot of the sgRNA representation")
qcpc <- qcpc + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size = 0.8), axis.text=element_text(size=15), axis.title=element_text(size=20), axis.line = element_line(size = 0),  legend.title=element_text(size=0), legend.text=element_text(size=20), plot.title = element_text(size=22))
qcpc <- qcpc + labs(x = "PC3 (13.7%)", y = "PC4 (10.6%)") 
qcpc <- qcpc + scale_shape_manual(values=c(21, 2, 5)) #+ ggtitle("PCA Su CD4 dataset, 3-mers")
pdf("pca_full_species_set2.pdf", height = 5, width = 13)
qcpc
dev.off()

# intersect 
pc$label <- as.character(metadata$class_intersect[match(pc$sample, metadata$`GenBank Accession`)])
qcpc <- ggplot(data=pc, aes(x=PC1, y=PC2)) + geom_point( aes(color=label),  size = 5,pch=1, stroke = 0.4) # + scale_fill_manual(values=c("#9399ff", "#a9fffd", "#46cdcf",  "#f67280", "#c06c84", "#6c5b7b", "#355c7d"))# + #+ ggtitle("Boxplot of the sgRNA representation")
qcpc <- qcpc + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size = 0.8), axis.text=element_text(size=15), axis.title=element_text(size=20), axis.line = element_line(size = 0),  legend.title=element_text(size=0), legend.text=element_text(size=20), plot.title = element_text(size=22))
qcpc <- qcpc + labs(x = "PC1 (24%)", y = "PC2 (19%)") 
qcpc <- qcpc + scale_shape_manual(values=c(21, 2, 5)) #+ ggtitle("PCA Su CD4 dataset, 3-mers")
pdf("pca_full_intersect.pdf", height = 5, width = 7)
qcpc
dev.off()

# pathogenic 
pc$label <- as.character(metadata$class_pathogenic[match(pc$sample, metadata$`GenBank Accession`)])

qcpc <- ggplot(data=pc, aes(x=PC1, y=PC2)) + geom_point( aes(color=label),  size = 5,pch=1, stroke = 0.4) # + scale_fill_manual(values=c("#9399ff", "#a9fffd", "#46cdcf",  "#f67280", "#c06c84", "#6c5b7b", "#355c7d"))# + #+ ggtitle("Boxplot of the sgRNA representation")
qcpc <- qcpc + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size = 0.8), axis.text=element_text(size=15), axis.title=element_text(size=20), axis.line = element_line(size = 0),  legend.title=element_text(size=0), legend.text=element_text(size=20), plot.title = element_text(size=22))
qcpc <- qcpc + labs(x = "PC1 (24%)", y = "PC2 (19%)") 
qcpc <- qcpc + scale_shape_manual(values=c(21, 2, 5)) #+ ggtitle("PCA Su CD4 dataset, 3-mers")
pdf("pca_full_class_pathogenic_set1.pdf", height = 5, width = 6)
qcpc
dev.off()

pc$label <- as.character(metadata$class_pathogenic[match(pc$sample, metadata$`GenBank Accession`)])
qcpc <- ggplot(data=pc, aes(x=PC3, y=PC4)) + geom_point( aes(color=label),  size = 5,pch=1, stroke = 0.4) # + scale_fill_manual(values=c("#9399ff", "#a9fffd", "#46cdcf",  "#f67280", "#c06c84", "#6c5b7b", "#355c7d"))# + #+ ggtitle("Boxplot of the sgRNA representation")
qcpc <- qcpc + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size = 0.8), axis.text=element_text(size=15), axis.title=element_text(size=20), axis.line = element_line(size = 0),  legend.title=element_text(size=0), legend.text=element_text(size=20), plot.title = element_text(size=22))
qcpc <- qcpc + labs(x = "PC3 (13.7%)", y = "PC4 (10.6%)") 
qcpc <- qcpc + scale_shape_manual(values=c(21, 2, 5)) #+ ggtitle("PCA Su CD4 dataset, 3-mers")
pdf("pca_full_class_pathogenic_set2.pdf", height = 5, width = 6)
qcpc
dev.off()
