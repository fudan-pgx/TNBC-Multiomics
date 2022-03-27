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
## cdefg. the consistency of variant numbers
FUSCCTNBC_var_num <- read.csv("./results/tables/FUSCCTNBC_var_num.csv",row.names = 1)
FUSCCTNBC_var_num.forplot <- t(FUSCCTNBC_var_num) %>% as.data.frame()
FUSCCTNBC_var_num.forplot$group <- rownames(FUSCCTNBC_var_num.forplot)
FUSCCTNBC_var_num.forplot$compare <- c("new","old","new","old","new","old","new","old","new","old")
FUSCCTNBC_var_num.forplot$dataset <- sapply(strsplit(as.character(FUSCCTNBC_var_num.forplot$group),"_"),
                                            function(x){x[1]})
var_num <- melt(FUSCCTNBC_var_num.forplot, id=c("dataset","compare"))
var_num.forplot <- dcast(var_num,variable+dataset~compare)

## outliers remove
outlier(FUSCCTNBC_var_num);outlier(FUSCCTNBC_var_num,opposite=TRUE)

### Remove the outliers in each column of the data, just enter the data (modify the name of dat), 
### you need to make sure that each column of the entered data needs to remove the extreme values
#dat is the input data, i.e. the data frame or matrix from which the extreme values need to be removed
## define the outlier remove function
outliers_remove <- function(Dat) {
  for (i in 1:ncol(Dat)) {
    Dat[, i] <- as.numeric(as.character(Dat[, i])) #Modify the target column to a numeric variable
    #extreme value upper bound : Q3+k(Q3-Q1)
    outlier_limup <- 3 * IQR(Dat[, i], na.rm = TRUE)+ quantile(Dat[, i], 3 / 4,na.rm = TRUE, names = FALSE)
    #3* Lower bound for extreme values : Q1-k(Q3-Q1)
    outlier_limdown <- quantile(Dat[, i], 1 / 4, na.rm = TRUE, names = FALSE)- 3 * IQR(Dat [, i], na.rm = TRUE)
    # Directly remove the row where the extreme value is located
    Dat <- Dat[!Dat[, i] >= outlier_limup & ! Dat[, i] <= outlier_limdown, ]
    # Change the extreme value point to null
    #Dat[Dat[, i] >= outlier_limup | Dat[, i] <= outlier_limdown, i] = ""
  }  
  return (Dat)
}

### single scatter plot
## VarScan
p1 <- ggscatter(outliers_remove(FUSCCTNBC_var_num[,c("VarScan_new","VarScan_old")]),
                x = "VarScan_new", y = "VarScan_old",
                add = "reg.line", conf.int = TRUE, color=mypal[3], alpha=0.5,
                add.params = list(linetype=2, color = "grey40", fill = "lightgray"))+
  stat_cor(method = "pearson",alternative = "two.sided",r.digits = 4,
           p.digits = 4,label.x = 100, label.y = 20)+
  xlim(0,500)+ylim(0,500)+
  labs(title="VarScan", x="", y="")+
  theme(plot.title = element_text(hjust = 0.5));p1

## TNscope
p2 <- ggscatter(outliers_remove(FUSCCTNBC_var_num[,c("TNscope_new","TNscope_old")]),
                x = "TNscope_new", y = "TNscope_old",
                add = "reg.line", conf.int = TRUE,color=mypal[4], alpha=0.5,
                add.params = list(linetype=2, color = "grey40", fill = "lightgray"))+
  stat_cor(method = "pearson",alternative = "two.sided",r.digits = 4,
           p.digits = 4,label.x = 250, label.y = 100)+
  xlim(0,500)+ylim(0,500)+
  labs(title="TNscope", x="", y="")+
  theme(plot.title = element_text(hjust = 0.5));p2

## TNseq
p3 <- ggscatter(outliers_remove(FUSCCTNBC_var_num[,c("TNseq_new","TNseq_old")]),
                x = "TNseq_new", y = "TNseq_old",
                add = "reg.line", conf.int = TRUE, color=mypal[5], alpha=0.5,
                add.params = list(linetype=2, color = "grey40", fill = "lightgray"))+
  stat_cor(method = "pearson",alternative = "two.sided",r.digits = 4,
           p.digits = 4,label.x = 100, label.y = 20)+
  xlim(0,500)+ylim(0,500)+
  labs(title="TNseq", x="", y="")+
  theme(plot.title = element_text(hjust = 0.5));p3

## Filter Mutation
p4 <- ggscatter(outliers_remove(FUSCCTNBC_var_num[,c("mutation_set_new","mutation_set_old")]),
                x = "mutation_set_new", y = "mutation_set_old",
                add = "reg.line", conf.int = TRUE, color=mypal[2], alpha=0.5,
                add.params = list(linetype=2, color = "grey40", fill = "lightgray"))+
  stat_cor(method = "pearson",alternative = "two.sided",r.digits = 4,
           p.digits = 4,label.x = 100, label.y = 20)+
  xlim(0,500)+ylim(0,500)+
  labs(title="Filter Mutation", x="", y="")+
  theme(plot.title = element_text(hjust = 0.5));p4

## Final Mutation
p5 <- ggscatter(outliers_remove(FUSCCTNBC_var_num[,c("final_new","final_old")]),
                x = "final_new", y = "final_old",
                add = "reg.line", conf.int = TRUE, color=mypal[1], alpha=0.5,
                add.params = list(linetype=2, color = "lightblue", fill = "lightgray"))+
  stat_cor(method = "spearman",alternative = "two.sided",r.digits = 4,
           p.digits = 4,label.x = 100, label.y = 20)+
  xlim(0,400)+ylim(0,400)+
  labs(title="Final Mutation", x="", y="")+
  theme(plot.title = element_text(hjust = 0.5));p5

p6 <- p5+p4+p1+p2+p3;p6

panel.cdefg <- ggarrange(p5,p4,p1,p2,p3,nrow = 1, labels = c("c","d","e","f","g")) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.cdefg
ggsave("./results/charts/panel cdefg.pdf",panel.cdefg,width=20,height=4,dpi=300)

## h. the consistency of variant genes
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

panel.h <- ggarrange(FUSCCTNBC_JI_plo, nrow = 1, labels = c("h")) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.h
ggsave("./results/charts/panel h.pdf",panel.h,width=6,height=4,dpi=300)

## i. the consistency of hotspot mutated genes
compare_genes_info <- read.csv("./results/tables/compare_genes_info.csv",row.names = 1)
compare_genes_info$percentile.x <- log10(compare_genes_info$percentile.x)
compare_genes_info$percentile.y <- log10(compare_genes_info$percentile.y)
# compare_genes_forplo <- compare_genes_info[-18,]

p7 <- ggscatter(compare_genes_info, x = "percentile.y", y = "percentile.x",
                add = "reg.line", conf.int = TRUE, color = "Var1",palette = pal_20,
                add.params = list(linetype=2, color = "lightblue", fill = "lightgray"))+
  stat_cor(method = "pearson",alternative = "two.sided",r.digits = 4,
           p.digits = 4,label.x = -1.2, label.y = -1.5)+
  labs(title="", x="Log10(Frequency) (New pipeline)", y="Log10(Frequency) (Old pipeline)")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")+mytheme+
  geom_text_repel(aes(percentile.x,percentile.y,
                       color=factor(Var1),label=Var1),
                   point.padding = NA, 
                   arrow = arrow(length=unit(0.01, "npc")),
                   max.overlaps = 20);p7

ggsave('./results/charts/hotspot_genes_precentile_scatter new.png',p7,width=5,height=5)

panel.i <- ggarrange(p7, nrow = 1, labels = c("i")) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.i
ggsave("./results/charts/panel i.pdf",panel.i,width=4,height=4,dpi=300)

panel.hi <- ggarrange(panel.h,panel.i,widths = c(1,1)) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.hi
ggsave("./results/charts/panel hi.pdf",panel.hi,width=12,height=4,dpi=300)

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

panel.j <- ggarrange(p8,p9, nrow = 1, labels = c("j"),widths = c(1,1)) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.j
ggsave("./results/charts/panel.j.pdf",panel.j,width=12,height=6,dpi=300)

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


panel.k <- ggarrange(FUSCCTNBC_CNV_JI_plo, nrow = 1, labels = c("k"),widths = c(2,4)) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.k
ggsave("./results/charts/panel.k.png",panel.k,width=8,height=6,dpi=300)
ggsave("./results/charts/panel.k.pdf",panel.k,width=8,height=6,dpi=300)

panel.jk <- ggarrange(panel.j, panel.k,nrow = 1, widths = c(2,3)) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.jk

fig4 <- ggarrange(panel.ab,panel.cdefg,panel.hi,panel.jk,nrow = 4) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));fig4
ggsave("./results/charts/fig4_0325.pdf",fig4,width=16,height=20,dpi=300)

## remove all
rm(list=ls())
