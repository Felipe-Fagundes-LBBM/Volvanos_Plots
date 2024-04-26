library(ggplot2)
library(ggrepel)

setwd('C:/Users/User/Desktop/adaptação_fontes_andre')

### dataset
de <- read.delim("volvano_andre/DESeq2_All_DEGs_SAMS_x_NonSAMS_ATV10uM_vs_VHC_comma.txt", sep= "\t", header = T, dec = ",")

de$diffexpressed <- "NO"
de$diffexpressed[de$log2FoldChange > 1 & de$pvalue < 0.00005] <- "UP"
de$diffexpressed[de$log2FoldChange < -1 & de$pvalue < 0.00005] <- "DOWN"
de$diffexpressed <- factor(de$diffexpressed, 
                          levels = c("UP","DOWN","NO"))
de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- de$gene_name[de$diffexpressed != "NO"]

### Plot

mycolors <- c("blue2","coral2","gray90")
names(mycolors) <- c("DOWN", "UP", "NO")

pdf("volvano_andre/Volcano_DEGs_SAM_nonSAM_ATV1.pdf")

ggplot(data = de , aes(x=log2FoldChange, y=-log10(pvalue))) +
  labs(title="", 
       x="log2 (Fold-Change)", 
       y="-log10 (P-value)")+
  scale_fill_manual(values=mycolors)+
  geom_point(aes(fill=diffexpressed), shape=21, size = 5, alpha = .9) +
  scale_x_continuous(breaks = c(-1,-3,-6,-9,1,3,6,9), limits = c(-9,9))+
  scale_y_continuous(breaks = c(0:9), limits = c(0,9))+
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.00005), col="red")+
  theme(axis.text = element_text(size = 24, color = "#636363"),
        axis.title = element_text(size = 30, color = "#636363"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        plot.title = element_text(size = 17, hjust = 0.5))+
  geom_label_repel(data= de, aes(label=delabel), color="black",
                   size=5, fontface = "bold", force = 8 )
dev.off()

###########        ROS

### dataset
de <- read.delim("volvano_andre/DESeq2_All_DEGs_SAMS_x_NonSAMS_ROS10uM_vs_VHC_comma (1).txt", sep= "\t", header = T, dec = ",")

de$diffexpressed <- "NO"
de$diffexpressed[de$log2FoldChange > 1 & de$pvalue < 0.00005] <- "UP"
de$diffexpressed[de$log2FoldChange < -1 & de$pvalue < 0.00005] <- "DOWN"
de$diffexpressed <- factor(de$diffexpressed, 
                           levels = c("UP","DOWN","NO"))
de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- de$gene_name[de$diffexpressed != "NO"]

### Plot

mycolors <- c("blue2","coral2","gray90")
names(mycolors) <- c("DOWN", "UP", "NO")

pdf("volvano_andre/Volcano_DEGs_SAM_nonSAM_ROS1.pdf")
ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue))) +
  labs(title="", 
       x="log2 (Fold-Change)", 
       y="-log10 (p-value)")+
  scale_fill_manual(values=mycolors)+
  geom_point(aes(fill=diffexpressed), shape=21, size = 5, alpha = .9) +
  scale_x_continuous(breaks = c(-1,-3,-6,-9,1,3,6,9), limits = c(-9,9))+
  scale_y_continuous(breaks = c(0:9), limits = c(0,9))+
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.00005), col="red")+
  theme(axis.text = element_text(size = 24, color = "#636363"),
        axis.title = element_text(size = 30, color = "#636363"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        plot.title = element_text(size = 17, hjust = 0.5))+
  geom_label_repel(data= de, aes(label=delabel), color="black",
                   size=4, fontface = "bold", force = 8 )
dev.off()


############proteomics ATV

limma_DEPs_All_ATV_SxVHC_S_x_ATV_NSxVHC_NS


### dataset
de <- read.delim("limma_DEPs_All_ATV_SxVHC_S_x_ATV_NSxVHC_NS.txt", sep= "\t", header = T, dec = ",")

de$diffexpressed <- "NO"
de$diffexpressed[de$logFC > 1 & de$PValue < 0.00005] <- "UP"
de$diffexpressed[de$logFC < -1 & de$PValue < 0.00005] <- "DOWN"
de$diffexpressed <- factor(de$diffexpressed, 
                           levels = c("UP","DOWN","NO"))
de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] = de$gene_symbol[de$diffexpressed!="NO"]

### Plot

mycolors <- c("blue2","coral2","gray90")
names(mycolors) <- c("DOWN", "UP", "NO")

pdf("Proteomics_Volcano_DEGs_SAM_nonSAM_ATV.pdf")
ggplot(data=de, aes(x=logFC, y=-log10(PValue))) +
  labs(title="", 
       x="log2 (Fold-Change)", 
       y="-log10 (p-value)")+
  scale_fill_manual(values=mycolors)+
  geom_point(aes(fill=diffexpressed), shape=21, size = 5, alpha = .9) +
  scale_x_continuous(breaks = c(-1,-3,-6,-9,1,3,6,9), limits = c(-9,9))+
  scale_y_continuous(breaks = c(0:9), limits = c(0,9))+
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.00005), col="red")+
  theme(axis.text = element_text(size = 20, color = "black"),
        axis.title = element_text(size = 26, color = "black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        plot.title = element_text(size = 17, hjust = 0.5))+
  geom_label_repel(data= de, aes(label=delabel), color="black",
                  size=4, fontface = "bold", force = 8 )
dev.off()

############proteomics ROS

### dataset
deros <- read.delim("limma_DEPs_All_ROS_SxVHC_S_x_ROS_NSxVHC_NS.txt", sep= "\t", header = T, dec = ",")

deros$diffexpressed <- "NO"
deros$diffexpressed[deros$logFC > 1 & deros$PValue < 0.00005] <- "UP"
deros$diffexpressed[deros$logFC < -1 & deros$PValue < 0.00005] <- "DOWN"
deros$diffexpressed <- factor(deros$diffexpressed, 
                           levels = c("UP","DOWN","NO"))
deros$delabel <- NA
deros$delabel[deros$diffexpressed != "NO"] = deros$gene_symbol[deros$diffexpressed!="NO"]

### Plot

mycolors <- c("blue2","coral2","gray90")
names(mycolors) <- c("DOWN", "UP", "NO")

pdf("Proteomics_Volcano_DEGs_SAM_nonSAM_ROS.pdf")
p=ggplot(data=deros, aes(x=logFC, y=-log10(PValue))) +
  labs(title="", 
       x="log2 (Fold-Change)", 
       y="-log10 (p-value)")+
  scale_fill_manual(values=mycolors)+
  geom_point(aes(fill=diffexpressed), shape=21, size = 5, alpha = .9) +
  scale_x_continuous(breaks = c(-1,-3,-6,-9,1,3,6,9), limits = c(-9,9))+
  scale_y_continuous(breaks = c(0:9), limits = c(0,9))+
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.00005), col="red")+
  theme(axis.text = element_text(size = 20, color = "black"),
        axis.title = element_text(size = 26, color = "black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        plot.title = element_text(size = 17, hjust = 0.5))+
  geom_label_repel(data= deros, aes(label=delabel), color="black",
                   size=4, fontface = "bold", force = 8 )
dev.off()
p


