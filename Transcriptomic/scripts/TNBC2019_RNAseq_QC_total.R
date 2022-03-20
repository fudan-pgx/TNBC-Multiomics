## TNBC2019 RNAseq QC
## Author：Qingwang Chen; Yaqing Liu; Ying Yu；Zhihui Li
## Date：2022-3-1
## Version：3.0

# Set work directory and import libraries
rm(list = ls (all = TRUE))

source("./scripts/libraries.R")
source("./scripts/parameter.R")

######### metadata import
qc_metadata <- read.csv("./data/TNBC2019_RNASeqQC_dir_metadata.csv",header=T)

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
# ggsave("./results/charts/Mapping_atioForSamples_RNAseq.png",width=8,height=6,dpi=300)

min(qc_MR$percentage_aligned) # 96.24614

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
######### 3. genomic region
genomic <- qc_metadata[,c("genomic_origin_id", "Exonic", "Intronic", 
                          "Intergenic")] %>% distinct(genomic_origin_id, .keep_all = T)

genomic$genomic_origin_id <- gsub(".percent","",genomic_region$genomic_origin_id,fix=T)
rownames(genomic) <- genomic$genomic_origin_id
genomic$genomic_origin_id <- NULL

genomic.forplot<-melt(as.matrix(genomic))
colnames(genomic.forplot)<-c("file","region","value")
genomic.forplot$region <- factor(genomic.forplot$region, levels=rev(levels(genomic.forplot$region)))

G_barplo<- ggplot(genomic.forplot,aes(x=file,y=value/1000,fill=region ))+
  geom_bar(position = "fill",stat="identity") +
  theme_bw()+
  mytheme+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom")+
  scale_fill_brewer(palette = "Paired")+
  labs(fill = "Region",x="",y="% Region percent")+
  scale_y_continuous(expand = c(0,0.01)); G_barplo

genomic$ExonicRatio<-apply(genomic[,1:3],1,function(x){round(x[1]*100/sum(x),2)})
genomic$IntronicRatio<-apply(genomic[,1:3],1,function(x){round(x[2]*100/sum(x),2)})
genomic$IntergenicRatio<-apply(genomic[,1:3],1,function(x){round(x[3]*100/sum(x),2)})

G_barplo_order <- G_barplo+scale_x_discrete(limits=rownames(genomic)[order(genomic$ExonicRatio)])
G_barplo_order

ggsave("./results/charts/gene_regionSamples_RNAseq.png",width=6 ,height=4,dpi=300)

panel.c <- ggarrange(G_barplo_order, nrow = 1,labels = c("c"))+
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.c

panel.bc <- ggarrange(panel.b,panel.c,widths = c(3,2))+
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.bc
ggsave("./results/charts/panel bc.pdf",panel.bc,width=12,height=4,dpi=300)


#------------------------------------------------------------------------------#
######### 4. Coverage
# subset the coverage part of QC metadata
qc_coverage <- qc_metadata[,c("sample_id", "mean_coverage",
                              "tissue_type")]  %>% distinct(sample_id, .keep_all = T)

## plot
Cov_histogram <- gghistogram(qc_coverage, x = "mean_coverage",
                             add = "median", rug = TRUE, 
                             color = "tissue_type", fill = "tissue_type", alpha = 0.5,
                             palette = mypal[c(5,4)], xlab = "Coverage", ylab = "Sample counts", 
                             legend.title = "Type") + mytheme + 
  theme(legend.position = "none", legend.margin=unit(-.05,"cm")); Cov_histogram

## export
ggsave("./results/charts/coverageForSamples_RNAseq.png",width=4,height=4,dpi=300)

aggregate(qc_coverage$mean_coverage, by = list(qc_coverage$tissue_type), FUN = mean)
# rm(list = ls(pattern = "coverage|histogram"))

#------------------------------------------------------------------------------#
######## 5. Insert size
qc_ins_size <- qc_metadata[,c("sample_id", "median_insert_size", 
                              "tissue_type")] %>% distinct(sample_id, .keep_all = T)
## plot
Ins_denplo_area <- ggdensity(qc_ins_size, x = "median_insert_size",
                             add = "median", rug = TRUE, 
                             color = "tissue_type", fill = "tissue_type", alpha = 0.5, 
                             palette = mypal[c(5,4)],xlab = "Insert Size", ylab = "Density", 
                             legend.title = "Type")+ mytheme + 
                   theme(legend.position = "none");Ins_denplo_area

ggsave("./results/charts/Insert_sizeForSamples_RNAseq.png",width=6,height=4,dpi=300)

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

ggsave("./results/charts/CG_content_vio_ForSamples_RNAseq.png",width=6,height=4,dpi=300)
aggregate(qc_gc_content$avg_gc, by = list(qc_gc_content$tissue_type), FUN = mean)

panel.def <- ggarrange(Cov_histogram,Ins_denplo_area, GC_vioplo, nrow = 1, widths = c(1,1,1), 
                     labels = c("d","e","f"))+theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.def
ggsave("./results/charts/panel def.pdf",panel.def,width=16,height=6,dpi=300)

#------------------------------------------------------------------------------#
## TNBC2019 RNAseq exp files QC
######## Data Clean
expr.mat.fpkm <- fread("./data/tnbc2019_448samples_dir_fpkm.csv") %>% as.data.frame()
rownames(expr.mat.fpkm) <- expr.mat.fpkm$gene_id
expr.mat.fpkm$gene_id <- NULL
expr.library <- colnames(expr.mat.fpkm) %>% as.data.table()
colnames(expr.library) <- "library"
expr.library$type<- sapply(strsplit(as.character(expr.library$library),"_| "),function(x){x[2]})
expr.library$type<-gsub("rep","TT",expr.library$type,fix=T)
expr.library$type[which(is.na(expr.library$type))] <- "TT"
expr.library$base_name <- sapply(strsplit(as.character(expr.library$library),"_| "),function(x){x[1]})
expr.library <- unite(expr.library, "colnames", base_name, type, sep = "_", remove = FALSE)
## Check if there are duplicates
sum(duplicated(expr.library$colnames))
colnames(expr.mat.fpkm) <- expr.library$colnames
saveRDS(expr.mat.fpkm,"./RData/expr_mat_fpkm_TNBC.rds")

## Add gene Name
gene_name <- read.table("./data/D5_1.gene.abundance.txt",header = T,row.names = NULL,sep = "\t")
expr.mat.fpkm$Gene_Name <- gene_name$Gene.Name[match(rownames(expr.mat.fpkm),gene_name$Gene.ID)]

## Calculate the mean value of this expression profile based on the gene name of expr.mat.fpkm
## If there are genes with the same name, take the max value of the identical genes
expr.mat.fpkm_g <- apply(expr.mat.fpkm[,1:448],2,function(x){tapply(x,expr.mat.fpkm$Gene_Name,max)})
dim(expr.mat.fpkm_g)
### [1] 56858   448
## Export exp files with gene name
saveRDS(expr.mat.fpkm_g,"./RData/gene_fpkm_r56858c448_genesymbol_TNBC2019_20220307.rds")

## Z-Score
expr.mat.fpkm_g_zscore <-  t(scale(t(expr.mat.fpkm_g)))
## Export Z-score exp files with gene name
saveRDS(expr.mat.fpkm_g_zscore,"./RData/Z_score_gene_fpkm_r56858c448_genesymbol_TNBC2019_20220307.rds")

#------------------------------------------------------------------------------#
######## Exp Data Process
## 1. Gender-checked sex cluster
### Import of gender-checked related genes and exp files
sexgene <- read.table("./data/sexgenelist.txt",header=T,sep="\t",stringsAsFactors = FALSE)
gene_fpkm <- readRDS("./Rdata/gene_fpkm_r56858c448_genesymbol_TNBC2019_20210819.rds")
dim(gene_fpkm)
### Calculate log(expr+1)
log_gene_expr <- apply(gene_fpkm,2,function(x){log2(x+1)})

### Screening sample for lines associated with sex-checked genes
logexpr_sex <- log_gene_expr[rownames(log_gene_expr) %in% sexgene$GeneSymbol,]

# plot
sex_heatmap <- pheatmap(logexpr_sex,scale = "column",distance_rows = "euclidean",
                        clustering_distance_cols = "euclidean",
                        colorRampPalette(c('#436eee','white','#EE0000'))(100),
                        clustering_method = "ward.D",show_rownames = TRUE,show_colnames=FALSE);sex_heatmap

# pheatmap to ggplot
sex_heatmap_g = as.ggplot(sex_heatmap)

ggsave("./results/charts/sex check total.pdf",sex_heatmap_g,width=6,height=4,dpi=300)

### Filter out samples with abnormal expression of RPS4Y1 gene
sample_gene_sex <- t(logexpr_sex) %>% as.data.frame()
max(sample_gene_sex$RPS4Y1)
# [1] 1.126515 FPKM.FUSCCTNBC475 (TT)
# rm(list=ls())

#------------------------------------------------------------------------------#
## 2. PCA clustering calibration
### Set group info 
group <- function(x){
  type<-rep("QC",length(x));
  type[grepl("PT",x)] <- "normal" ;
  type[grepl("TT",x)] <- "tumor";
  type <- factor(type,ordered=TRUE,levels=c("tumor","normal"))
  type
}

expr.mat.fpkm<-readRDS("./RData/expr_mat_fpkm_TNBC.rds")
# expr.mat.fpkm <- readRDS("./Rdata/gene_fpkm_r56858c448_genesymbol_TNBC2019_20210819.rds")

group.info <- data.table(colnames(expr.mat.fpkm),group(colnames(expr.mat.fpkm)))
colnames(group.info) <- c("sample_ID","group")

### 30% filter 
filtered.expr.mat.ft <- expr.mat.fpkm[apply(expr.mat.fpkm, MARGIN=1, 
                                            FUN=function(x) {(sum(x==0)/length(x))<=0.3}),]
log_filtered.expr.mat.ft <- apply(filtered.expr.mat.ft,2,function(x){log2(x+1)})
dim(log_filtered.expr.mat.ft)
## Filtered data PCA analysis
pca.all <- prcomp(t(log_filtered.expr.mat.ft),retx=T)

### Create horizontal and vertical labels
# PC3comp<-summary(pca.all)$importance[2,3]*100

group.info$group <- factor(group.info$group,levels = c("normal","tumor"))
pca_cluster <- fviz_pca_ind(pca.all, 
                            label = "none",
                            addEllipses = T,
                            geom.ind = c("point", "text"),
                            repel = T,
                            habillage = group.info$group,
                            ellipse.level=0.95,
                            palette = mypal[c(5,4)])+
  theme(axis.text = element_text(size=16,face='bold',color='gray40'),
        axis.title = element_text(size=16,face='bold'),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16,color='gray40'),
        plot.subtitle = element_text(hjust=0.5,size=16),
        plot.title = element_text(hjust=0.5,size=16))+
  scale_shape_manual(values = c(15:16))+
  labs(x = sprintf("PC1 (%.2f%%)", summary(pca.all)$importance[2,1]*100),
       y = sprintf("PC2 (%.2f%%)", summary(pca.all)$importance[2,2]*100),
       title=paste('PCA-All ',ncol(log_filtered.expr.mat.ft),'samples',' N=',nrow(log_filtered.expr.mat.ft)))+
  theme(aspect.ratio = 1/1);pca_cluster

ggsave("./results/charts/PCA_all_448samples.pdf",pca_cluster,width=6,height=4,dpi=300)

panel.gh <- ggarrange(sex_heatmap_g, pca_cluster, nrow = 1, widths = c(1, 1), 
                     labels = c("e","f"))+theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.gh
ggsave("./results/charts/panel gh.pdf",panel.gh,width=10,height=4,dpi=300)


### FASTQC
img <- readPNG("./data/fastqc.png")
# plot with picture as layer
fastqc_p <- ggplot(mapping = aes(1:12, 1:12)) +
  annotation_raster(img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
  theme_bw() + xlab("")+ylab("")+
  theme(panel.border = element_blank());fastqc_p

panel.a <- ggarrange(fastqc_p,nrow = 1,labels = c("a"))+
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm")); panel.a

fig2 <- ggarrange(panel.a,panel.bc,panel.def,panel.gh,nrow = 4,heights = c(1,1,1,1))+
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));fig2
ggsave("./results/charts/fig2.pdf",fig2,width=12,height=16.5,dpi=300)


######## Delete unrelated variables
rm(list = ls())
