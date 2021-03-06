## Title: TNBC Bridge Figures Merge  ##
## Author: Qingwang Chen             ##
## Date: 2022-03-02                  ##  
## Version: V1.0                     ##
#######################################

# Set work directory and import libraries
rm(list = ls (all = TRUE))
source("./scripts/libraries.R")
source("./scripts/parameter.R")

#### --------------------------------RNAseq--------------------------------#### 
## b. FC Correlation
compare_FC <- read.csv("./results/tables/compare_FC.csv")

compare_FC_plot <- ggscatter(compare_FC, x = "log2FC.x", y = "log2FC.y",
                             add = "reg.line", conf.int = TRUE, alpha = 0.1,
                             add.params = list(linetype=2, color = "lightblue", fill = "lightgray"))+
  stat_cor(method = "pearson",alternative = "two.sided",r.digits = 4,
           p.digits = 4,label.x = 0, label.y = -2.5)+
  labs(title="", x="Log2FC (New)", y="Log2FC (Old)")+
  theme(plot.title = element_text(hjust = 0.5));compare_FC_plot

panel.b <- ggarrange(compare_FC_plot, nrow = 1, labels = c("b")) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.b
ggsave("./results/charts/panel b.pdf",panel.b,width=6,height=6,dpi=300)

panel.ab <- ggarrange(NULL,panel.b,widths = c(2,1)) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.ab
ggsave("./results/charts/panel ab.pdf",panel.ab,width=12,height=4,dpi=300)

#### ----------------------------------WES----------------------------------#### 
## c. the consistency of variant genes
### whole genes
FUSCCTNBC_JI <- read.csv("./results/tables/FUSCCTNBC_JI.csv",row.names = 1)
FUSCCTNBC_JI_forplot <- melt(FUSCCTNBC_JI)
FUSCCTNBC_JI_forplot$variable <- factor(FUSCCTNBC_JI_forplot$variable, 
                                        levels = c("Final.Mutation","Filter.Mutation",
                                                   "VarScan","TNscope","TNseq"))

FUSCCTNBC_JI_plo <- ggplot(FUSCCTNBC_JI_forplot,aes(x=variable, y=value,
                                                    color=variable,fill=variable))+
  geom_jitter(alpha = 0.5,size = 0.8)+
  geom_boxplot(alpha = 0.5,width = 0.8,
               outlier.shape = NA)+
  scale_fill_npg()+
  scale_color_npg()+
  theme_bw()+mytheme+
  labs(title="", x="", y="Jaccard Index")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylim(0,1);FUSCCTNBC_JI_plo
ggsave('./results/charts/FUSCCTNBC_JI_boxplot.png',FUSCCTNBC_JI_plo,width=4,height=3)

panel.c <- ggarrange(FUSCCTNBC_JI_plo, nrow = 1, labels = c("c")) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.c
ggsave("./results/charts/panel c.pdf",panel.c,width=6,height=4,dpi=300)

## d. the consistency of hotspot mutated genes
compare_genes_info <- read.csv("./results/tables/compare_genes_info.csv",row.names = 1)
# compare_genes_forplo <- compare_genes_info[-18,]

p1 <- ggscatter(compare_genes_info, x = "percentile.y", y = "percentile.x",
                add = "reg.line", conf.int = TRUE, color = "Var1",palette = pal_20,
                add.params = list(linetype=2, color = "lightblue", fill = "lightgray"))+
  labs(title="", x="Frequency (New pipeline)", y="Frequency (Old pipeline)")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")+mytheme+
  geom_text_repel(aes(percentile.x,percentile.y,
                       color=factor(Var1),label=Var1),
                   point.padding = NA, 
                   arrow = arrow(length=unit(0.01, "npc")),
                   max.overlaps = 20);p1
# R = 0.993, p = 3.78e-16
p2 <- ggpar(p1,xscale = "log10",yscale = "log10");p2

ggsave('./results/charts/hotspot_genes_precentile_scatter new.png',p2,width=5,height=5)

panel.d <- ggarrange(p2, nrow = 1, labels = c("d")) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.d
ggsave("./results/charts/panel d.pdf",panel.d,width=4,height=4,dpi=300)

panel.cd <- ggarrange(panel.c,panel.d,widths = c(2,1)) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.cd
ggsave("./results/charts/panel cd.pdf",panel.cd,width=12,height=4,dpi=300)

#### ----------------------------------CNV----------------------------------#### 
# the consistency of CNV Peaks(venn.diagram)
amp_inter_venn <- read.csv("./data/amp_inter_venn.csv",header = T)
del_inter_venn <- read.csv("./data/del_inter_venn.csv",header = T)

data1 <- list()
data1$Amp_peak_new <- amp_inter_venn[,4]
data1$Amp_peak_old <- amp_inter_venn[,2]

data2 <- list()
data2$Del_peak_new <- del_inter_venn[,4]
data2$Del_peak_old <- del_inter_venn[,2]

venn.plot <- venn.diagram(
  x = list(Amp_peak_new = data1$Amp_peak_new,
           Amp_peak_old = data1$Amp_peak_old),
  filename = NULL,
  fill = c("#be482c", "#be1b32"),
  alpha = 0.6,
  label.col = "white",
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("#be482c", "#be1b32"),
  cat.cex = 2,
  cat.fontfamily = "serif",
  cat.fontface = "bold",
  margin = 0.05,
  cat.dist = c(0.03, 0.03),
  cat.pos = c(360, 180),
  rotation.degree=90
)
p8 <- as.ggplot(plot_grid(grobTree(venn.plot)));p8

venn.plot <- venn.diagram(
  x = list(Del_peak_new = data2$Del_peak_new,
           Del_peak_old = data2$Del_peak_old),
  filename = NULL,
  fill = c("#64999f", "#93b18d"),
  alpha = 0.6,
  label.col = "white",
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("#64999f", "#93b18d"),
  cat.cex = 2,
  cat.fontfamily = "serif",
  cat.fontface = "bold",
  margin = 0.05,
  cat.dist = c(0.03, 0.03),
  cat.pos = c(360, 180),
  rotation.degree=90
)
p9 <- as.ggplot(plot_grid(grobTree(venn.plot)));p9

panel.e <- ggarrange(p8,p9, nrow = 1, labels = c("e"),widths = c(1,1)) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.e
ggsave("./results/charts/panel.e.pdf",panel.e,width=12,height=6,dpi=300)

## the consistency of CNV genes (venn.diagram)
# new_gene_cnv <- read.csv("./results/tables/new_gene_cnv.csv",header = T,row.names = 1)
# old_gene_cnv <- read.csv("./results/tables/old_gene_cnv.csv",header = T,row.names = 1)

# data <- list()
# data$newgene <- rownames(new_gene_cnv)
# data$oldgene <- rownames(old_gene_cnv)
# 
# venn.plot <- venn.diagram(
#   x = list(new_gene = data$newgene,
#            old_gene = data$oldgene),
#   filename = NULL,
#   fill = c("#f39b7f", "#3c5488"),
#   alpha = 0.6,
#   label.col = "white",
#   cex = 1.5,
#   fontfamily = "serif",
#   fontface = "bold",
#   cat.col = c("#f39b7f", "#3c5488"),
#   cat.cex = 2,
#   cat.fontfamily = "serif",
#   cat.fontface = "bold",
#   margin = 0.05,
#   cat.dist = c(0.03, 0.03),
#   cat.pos = c(0, 180),
#   rotation.degree=270
# );grid.draw(venn.plot)
# p10 <- as.ggplot(plot_grid(grobTree(venn.plot)));p10
# ggsave("./results/charts/p10.pdf",p10,width=6,height=6,dpi=300)

## the consistency of CNV genes (Jaccard Index)
CNV_ji <- read.csv("./results/tables/CNV_ji.csv",header = T,row.names = 1,)
CNV_ji_forplot <-melt(CNV_ji)
CNV_ji_forplot$variable <- factor(CNV_ji_forplot$variable, 
                                        levels = c("Total.Genes","Netural","Gain",
                                                   "Loss"))
## Comparison of old and new processes
FUSCCTNBC_CNV_JI_plo <- ggplot(CNV_ji_forplot,aes(x=variable, y=value,
                                              color=variable,fill=variable))+
  # geom_jitter(alpha = 0.5,size = 1)+
  geom_boxplot(alpha = 0.5,width = 0.8,
               outlier.shape = NA)+
  scale_fill_aaas()+
  scale_color_aaas()+
  theme_bw()+mytheme+
  labs(title="", x="", y="Jaccard Index")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")+theme(axis.text.x = element_text(angle = 45, hjust = 1,size =12))+theme(axis.text.y = element_text(size =12))+
  ylim(0,1);FUSCCTNBC_CNV_JI_plo


panel.f <- ggarrange(FUSCCTNBC_CNV_JI_plo, nrow = 1, labels = c("f"),widths = c(2,4)) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.f
ggsave("./results/charts/panel.f.png",panel.f,width=8,height=6,dpi=300)
ggsave("./results/charts/panel.f.pdf",panel.f,width=8,height=6,dpi=300)

panel.ef <- ggarrange(panel.e, panel.f,nrow = 1, widths = c(2,3)) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.ef

fig4 <- ggarrange(panel.ab,panel.cd,panel.ef,nrow = 3) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));fig4
ggsave("./results/charts/fig4_0330.pdf",fig4,width=16,height=15,dpi=300)

## remove all
rm(list=ls())
