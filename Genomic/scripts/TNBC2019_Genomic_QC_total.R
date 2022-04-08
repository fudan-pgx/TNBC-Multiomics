## TNBC2019 WES QC
## Author：Qingwang Chen; Yaqing Liu
## Date：2021-8-11
## Verision：2.0

# 设置工作目录并导入相关库
rm(list = ls (all = TRUE))

source("./scripts/libraries.R")
source("./scripts/parameter.R")
source("./scripts/theme_color.R")

######### meatadata import
qc_metadata <- read.csv("./data/TNBC2019_WESQC_metadata.csv",header=T)

#------------------------------------------------------------------------------#
######### 1. Mapping ratio
qc_MR <- qc_metadata[,c("sample_id", "percentage_aligned",
                        "tissue_type")] %>% distinct(sample_id, .keep_all = T)
qc_MR$species <- "Human"

# MR_boxplo <- ggboxplot(qc_MR, x = "species", y = "percentage_aligned",
#                        add = c("median","jitter"),
#                        color = "species",
#                        palette = mypal[4],
#                        xlab = "", ylab = "% Mapping Ratio of Human Genome") + 
#   mytheme + theme(legend.position="none",axis.text.x = element_text(angle = 45, hjust = 1))+
#   ylim(0, 100); MR_boxplo
# ggsave("./results/charts/Mapping_atioForSamples_WES.png",width=8,height=6,dpi=300)

min(qc_MR$percentage_aligned) # 97.63026


#------------------------------------------------------------------------------#
######### 2. FastQ Screen (Contamination)
# colnames(qc_metadata)[grep("percentage",colnames(qc_metadata))]
fastqscreen <- qc_metadata[,c("fastqscreen_id", "Mouse",  "Yeast", "EColi", 
                              "Virus", "Vector")]

fastqscreen$fastqscreen_id <- gsub("_screen","",fastqscreen$fastqscreen_id,fix=T)
fastqscreen$sample_id <- sapply(strsplit(as.character(fastqscreen$fastqscreen_id),"_"),
                                function(x){paste(x[1],x[2],sep="_")})
fastqscreen$sample_id <- gsub("R1|R2|rep","TT",fastqscreen$sample_id)
fastqscreen[,2:6] <- fastqscreen[,2:6]/100
human_MR <- qc_MR[,c(1,2)]
MR <- merge(human_MR,fastqscreen,by="sample_id")  
rownames(MR)<-MR$fastqscreen_id 
MR$fastqscreen_id<-NULL; MR$sample_id <- NULL
colnames(MR)[1] <- "Human"

fastqscreen.forplot<-melt(as.matrix(MR))
colnames(fastqscreen.forplot)<-c("file","source","value")
fastqscreen.forplot$sample<-sapply(strsplit(as.character(fastqscreen.forplot$file),"_"),
                                   function(x){paste(x[1],x[2],sep="_")})
fastqscreen.forplot$type <- ifelse(fastqscreen.forplot$source == "Human","Homo sapiens","Other genomes")

###### plot
Con_boxplo <- ggboxplot(fastqscreen.forplot, x = "source", y = "value",
                        add = c("median","jitter"),
                        color = "source",
                        palette = mypal,
                        xlab = "", ylab = "% Mapping Ratio",
                        legend.title = "contamination") + mytheme+ 
  theme(legend.position="none",axis.text.x = element_text(angle = 45, hjust = 1))+
  ylim(0, 100)+facet_grid(.~type,scales = "free",space = "free"); Con_boxplo
ggsave("./results/charts/contaminForSamples_RNAseq.png",width=6,height=4,dpi=300)

panel.b <- ggarrange(Con_boxplo, nrow = 1, labels = c("b")) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.b
ggsave("./results/charts/panel b.pdf",panel.b,width=8,height=6,dpi=300)

#------------------------------------------------------------------------------#
######### 3. Coverage
# subset the coverage part of QC metadata
qc_coverage <- qc_metadata[,c("sample_id", "mean_coverage",
                              "tissue_type")]  %>% distinct(sample_id, .keep_all = T)

## plot
Cov_histogram <- gghistogram(qc_coverage, x = "mean_coverage",
                             add = "median", rug = TRUE, 
                             color = "tissue_type", fill = "tissue_type", alpha = 0.5,
                             palette = mypal[c(5,4)], xlab = "Coverage", ylab = "Sample counts", 
                             legend.title = "Type") + mytheme + 
  theme(legend.position = "bottom", legend.margin=unit(-.05,"cm")); Cov_histogram

## export
ggsave("./results/charts/coverageForSamples_WES.png",width=4,height=4,dpi=300)

aggregate(qc_coverage$mean_coverage, by = list(qc_coverage$tissue_type), FUN = mean)
# rm(list = ls(pattern = "coverage|histogram"))

panel.c <- ggarrange(Cov_histogram, nrow = 1, labels = c("c")) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.c
ggsave("./results/charts/panel c.pdf",panel.b,width=8,height=6,dpi=300)

panel.bc <- ggarrange(panel.b,panel.c, nrow = 1, widths = c(2, 1)) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.bc
ggsave("./results/charts/panel bc.pdf",panel.bc,width=12,height=4,dpi=300)

#------------------------------------------------------------------------------#
######## 4. Duplicate
qc_duplicate <- qc_metadata[,c("sample_id", "duplicate_picard", 
                               "tissue_type")] %>% distinct(sample_id, .keep_all = T)
## plot
qc_duplicate$tissue_type <- factor(qc_duplicate$tissue_type,levels = c("normal","tumor"))

Dup_vioplo <- ggviolin(qc_duplicate, "tissue_type", "duplicate_picard",
                       fill = "tissue_type", palette = mypal[c(5,4)], alpha = 0.5, 
                       add.params = list(fill = "white"),
                       add = "boxplot",
                       ylab = "% Duplication", xlab = "Tissue Type") + mytheme+
  theme(legend.position = "none"); Dup_vioplo

ggsave("./results/charts/duplicate_vio_ForSamples.png",width=4,height=4,dpi=300)

#------------------------------------------------------------------------------#
######## 5. Insert size
qc_ins_size <- qc_metadata[,c("sample_id", "median_insert_size", 
                              "tissue_type")] %>% distinct(sample_id, .keep_all = T)
## plot
Ins_denplo_area <- ggdensity(qc_ins_size, x = "median_insert_size",
                             add = "median", rug = TRUE, 
                             color = "tissue_type", fill = "tissue_type", alpha = 0.5, 
                             palette = mypal[c(5,4)],xlab = "Insert Size (bp)", ylab = "Density", 
                             legend.title = "Type")+ mytheme + 
  theme(legend.position = "none");Ins_denplo_area

ggsave("./results/charts/Insert_sizeForSamples_WES.png",width=6,height=4,dpi=300)

aggregate(qc_ins_size$median_insert_size, by = list(qc_ins_size$tissue_type), FUN = median)

#------------------------------------------------------------------------------#
######## 6. CG content
qc_gc_content <- qc_metadata[,c("sample_id", "avg_gc", 
                                "tissue_type")] %>% distinct(sample_id, .keep_all = T)
qc_gc_content$tissue_type <- factor(qc_gc_content$tissue_type,levels = c("normal","tumor"))

GC_vioplo <- ggviolin(qc_gc_content, "tissue_type", "avg_gc",
                      fill = "tissue_type",add.params = list(fill = "white"), alpha = 0.5,
                      palette = mypal[c(5,4)],
                      add = "boxplot",
                      ylab = "% GC", xlab = "Tissue Type")+ mytheme + 
  theme(legend.position = "none"); GC_vioplo

ggsave("./results/charts/CG_content_vio_ForSamples_WES.png",width=6,height=4,dpi=300)
aggregate(qc_gc_content$avg_gc, by = list(qc_gc_content$tissue_type), FUN = mean)

panel.def <- ggarrange(Ins_denplo_area,GC_vioplo,Dup_vioplo, nrow = 1, widths = c(1,1,1), 
                     labels = c("d","e","f"))+theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.def
ggsave("./results/charts/panel def.pdf",panel.def,width=12,height=6,dpi=300)

### FASTQC
img <- readPNG("./data/fastqc.png")
# plot with picture as layer
fastqc_p <- ggplot(mapping = aes(1:12, 1:12)) +
  annotation_raster(img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
  theme_bw() + xlab("")+ylab("")+
  theme(panel.border = element_blank());fastqc_p

panel.a <- ggarrange(fastqc_p,nrow = 1,labels = c("a"))+
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm")); panel.a




#------------------------------------------------------------------------------#
######## 7. OncoScan QC
OncoScan_CNV_QC <- read.csv("./data/OncoScan_CNV_QC220225.csv",header = T)
OncoScan_CNV_QC$MAPD_Con <- ifelse(OncoScan_CNV_QC$MAPD <= 0.3,"Pass","Fail")
# OncoScan_CNV_QC$ndWavinessSD_Con <- ifelse(OncoScan_CNV_QC$ndWavinessSD <= 0.3,"Pass","Fail")
OncoScan_CNV_QC$ndSNPQC_Con <- ifelse(OncoScan_CNV_QC$ndSNPQC >= 26,"Pass","Fail")

MAPD_histogram <- gghistogram(OncoScan_CNV_QC, x = "MAPD", rug = FALSE, 
            color = "MAPD_Con", fill = "MAPD_Con", alpha = 0.5,
            palette = mypal[c(8,7)], xlab = "MAPD Score", ylab = "Sample counts", 
            legend.title = "Type") + 
  mytheme + geom_vline(aes(xintercept=0.3), colour="#990000", linetype="dashed")+
  theme(legend.position = "none"); MAPD_histogram

# ndWavinessSD_histogram <- gghistogram(OncoScan_CNV_QC, x = "ndWavinessSD", rug = FALSE, 
#                                      color = "ndWavinessSD_Con", fill = "ndWavinessSD_Con", alpha = 0.5,
#                                      palette = mypal[c(8,7)], xlab = "ndWavinessSD Score", ylab = "", 
#                                      legend.title = "Type") + 
#   mytheme + geom_vline(aes(xintercept=0.3), colour="#990000", linetype="dashed")+
#   theme(legend.position = "none"); ndWavinessSD_histogram

ndSNPQC_histogram <- gghistogram(OncoScan_CNV_QC, x = "ndSNPQC", rug = FALSE, 
                                      color = "ndSNPQC_Con", fill = "ndSNPQC_Con", alpha = 0.5,
                                      palette = mypal[c(8,7)], xlab = "ndSNPQC Score", ylab = "Sample counts", 
                                      legend.title = "Type") + 
  mytheme + geom_vline(aes(xintercept=26), colour="#990000", linetype="dashed")+
  theme(legend.position = "none"); ndSNPQC_histogram

panel.gh <- ggarrange(MAPD_histogram,ndSNPQC_histogram,nrow = 1,labels = c("g","h"))+
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm")); panel.gh
ggsave("./results/charts/panel gh.pdf",panel.gh,width=6,height=4,dpi=300)

fig3 <- ggarrange(panel.a,panel.bc,panel.def,panel.gh,nrow = 4, widths = c(2,2,2,2))+
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));fig3
ggsave("./results/charts/fig3 0401.pdf",fig3,width=12,height=16,dpi=300)

######## Delete unrelated variables
rm(list = ls())

