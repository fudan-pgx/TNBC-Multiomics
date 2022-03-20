## Title: TNBC Bridge RNAseq Analysis##
## Author: Qingwang Chen             ##
## Date: 2022-03-02                  ##  
## Version: V2.0                     ##
#######################################

# Set work directory and import libraries
rm(list = ls (all = TRUE))
source("./scripts/libraries.R")
source("./scripts/parameter.R")

# ----------------------FC Correlation---------------------- #
# 1. Data process
load("./data/old_data_old_pipeline_rnaseq/FUSCCTNBC_RNAseqShi.Complete_log2_448x45308_V15_190209.Rdata")
fpkm_old_pipe <- FUSCCTNBC_RNAseqShi.Complete_log2

fpkm_new_pipe <- fread("./data/old_data_new_pipeline_rnaseq/tnbc2019_448samples_dir_fpkm.csv")
geneid2name <- fread("./data/geneid2name.csv")
## gene id to gene name
fpkm_new_pipe <- merge(geneid2name,fpkm_new_pipe)
fpkm_new_pipe$gene_id <- NULL
## max the fpkm of the same genes
fpkm_new_pipe <- aggregate(fpkm_new_pipe[,-1], list(fpkm_new_pipe$gene_symbol), FUN=max)
row.names(fpkm_new_pipe)<- fpkm_new_pipe$Group.1
fpkm_new_pipe$Group.1 <- NULL
colnames(fpkm_new_pipe) <- gsub("FPKM.","",colnames(fpkm_new_pipe))
## save the raw fpkm witn annotation
# saveRDS(fpkm_new_pipe,"./data/FUSCCTNBC_new_fpkm_448x56858_220302.rds")
## filer the genes whose fpkm=0 in 30% samples
filtered_fpkm_new_pipe <- fpkm_new_pipe[apply(fpkm_new_pipe, MARGIN=1, 
                                              FUN=function(x) {(sum(x==0)/length(x))<=0.3}),]
### 数据log2处理
filtered_fpkm_new_pipe_log2<-log2(filtered_fpkm_new_pipe+0.01)
dim(filtered_fpkm_new_pipe_log2)
## save the raw fpkm witn annotation
# saveRDS(filtered_fpkm_new_pipe_log2,"./data/FUSCCTNBC_new_filtered_fpkm_log2_448x25589_220302.rds")

# 3. FC & p-value calculation (88*2+272)
fpkm_new_pipe <- readRDS("./data/FUSCCTNBC_new_fpkm_448x56858_220302.rds")
filtered_fpkm_new_pipe <- fpkm_new_pipe[apply(fpkm_new_pipe, MARGIN=1, 
                                              FUN=function(x) {(sum(x==0)/length(x))<=0.3}),]

filtered_fpkm_new_pipe_PT <- filtered_fpkm_new_pipe[,grep("PT",colnames(filtered_fpkm_new_pipe))]
filtered_fpkm_new_pipe_TT <- filtered_fpkm_new_pipe[,-grep("PT",colnames(filtered_fpkm_new_pipe))]

fpkm_old_pipe <- 2^fpkm_old_pipe-1
filtered_fpkm_old_pipe <- fpkm_old_pipe[apply(fpkm_old_pipe, MARGIN=1, 
                                              FUN=function(x) {(sum(x==0)/length(x))<=0.3}),]
# filtered_fpkm_old_pipe <- unique(filtered_fpkm_old_pipe)
filtered_fpkm_old_pipe_PT <- filtered_fpkm_old_pipe[,grep("PT",colnames(filtered_fpkm_old_pipe))]
filtered_fpkm_old_pipe_TT <- filtered_fpkm_old_pipe[,-grep("PT",colnames(filtered_fpkm_old_pipe))]

## define FC function
foldChange<-function(exp_TT,exp_PT)
{
  exp_mean_TT <- apply(as.matrix(exp_TT),1,mean)
  exp_mean_PT <- apply(as.matrix(exp_PT),1,mean)
  exp_FC <- log(exp_mean_TT/exp_mean_PT,base=2)
  return(exp_FC)
}

new_pipe_FC <- foldChange(filtered_fpkm_new_pipe_TT,filtered_fpkm_new_pipe_PT)
new_pipe_FC_df <- as.data.frame(new_pipe_FC)
new_pipe_FC_df$gene_symbol <- rownames(new_pipe_FC_df)
old_pipe_FC <- foldChange(filtered_fpkm_old_pipe_TT,filtered_fpkm_old_pipe_PT)
old_pipe_FC_df <- as.data.frame(old_pipe_FC)
old_pipe_FC_df$gene_symbol <- rownames(old_pipe_FC_df)

## compare the FC of different pipelines
compare_FC <- merge(new_pipe_FC_df,old_pipe_FC_df,by="gene_symbol")
plot(compare_FC$new_pipe_FC,as.numeric(compare_FC$old_pipe_FC))  
#### cor=0.9512
cor.test(compare_FC$new_pipe_FC,as.numeric(compare_FC$old_pipe_FC))
write.csv(compare_FC,"./results/tables/compare_FC.csv")

## plot
compare_FC_plot <- ggscatter(compare_FC, x = "new_pipe_FC", y = "old_pipe_FC",
                             add = "reg.line", conf.int = TRUE, 
                             add.params = list(linetype=2, color = "lightblue", fill = "lightgray"))+
                   stat_cor(method = "pearson",alternative = "two.sided",r.digits = 4,
                            p.digits = 4,label.x = 2.5, label.y = -2.5)+
  labs(title="Compare of the log2FC of pipelines", x="New", y="Old")+
  theme(plot.title = element_text(hjust = 0.5));compare_FC_plot

ggsave('./results/charts/compare_FC_plot_scatter.png',compare_FC_plot,width=5,height=5)

## cal t-test
## define t-test function
T_test<-function(exp_total,exp_TT,exp_PT)
{
  ttest_result <- matrix(NA,dim(exp_total)[1],4)
  rownames(ttest_result) <- rownames(exp_total)
  for(i in (1:dim(exp_total)[1]))
  {
   ttest_result[i,4] <- shapiro.test(exp_total[i,])$p.value   #正态性检验  https://blog.csdn.net/weixin_34185512/article/details/92157419  
    temp <- t.test(exp_TT[i,],exp_PT[i,],alternative = "two.sided", pair = F,var.equal = T)
    ttest_result[i,1] <- temp$statistic 
    ttest_result[i,2] <- temp$p.value 
    print(i)
  }
  ttest_result[,3] <- p.adjust(ttest_result[,2],method = "BH",dim(exp_total)[1])
  colnames(ttest_result) <- c("statistic","p-value","adjusted-pvalue","shapiro")
  return(ttest_result)
}

### new_pipeline
ttest_result_new_pipe <- T_test(as.matrix(filtered_fpkm_new_pipe),
                       as.matrix(filtered_fpkm_new_pipe_TT),
                       as.matrix(filtered_fpkm_new_pipe_PT)) %>% as.data.frame()
ttest_result_new_pipe$log2FC <- new_pipe_FC

write.csv(ttest_result_new_pipe,"./results/tables/ttest_result_new_pipe.csv")

### old_pipeline
ttest_result_old_pipe <- T_test(as.matrix(filtered_fpkm_old_pipe),
                                as.matrix(filtered_fpkm_old_pipe_TT),
                                as.matrix(filtered_fpkm_old_pipe_PT)) %>% as.data.frame()
ttest_result_old_pipe$log2FC <- old_pipe_FC

write.csv(ttest_result_old_pipe,"./results/tables/ttest_result_old_pipe.csv")

# ----------------------Subtyping---------------------- #
## method from Dr. Ma (method before)
### filter the protein coding genes
Gene_annotation_TNBC <- read.csv("./data/RNA_seq_Gene_annotation_FUSCCTNBC.csv", header = T)
row.names(Gene_annotation_TNBC) <- Gene_annotation_TNBC$GeneSymbol
Gene_annotation_TNBC <- Gene_annotation_TNBC[which(!is.na(Gene_annotation_TNBC$Category)), ]
ProteinCoding <- Gene_annotation_TNBC$GeneSymbol[Gene_annotation_TNBC$Category == "protein_coding"]

### Expression dataset for subsequent test
# filtered_fpkm_old_pipe_TT 
logexpr_old_TT <- log2(filtered_fpkm_old_pipe_TT +1)
Tes_Expr <- logexpr_old_TT[ProteinCoding, ]
dim(Tes_Expr)

### Select genes for clustering (sd top 2000)
SDs <- NULL
for(i in 1:nrow(Tes_Expr))
{
  SDs <- c(SDs, sd(as.numeric(Tes_Expr[i,])))
  # if(i/1000 == round(i/1000)) {print(i/nrow(Tes_Expr))}
}
names(SDs) <- row.names(Tes_Expr)
SDs        <- sort(SDs, decreasing = T)
Clustering_Genes <- names(SDs[1:2000])
### 2000 genes

### Expression dataset for clustering
Tes_Expr <- Tes_Expr[Clustering_Genes, ]
dim(Tes_Expr)

### clustering by k-means
set.seed(1234)
km1 <- kmeans(t(Tes_Expr),4,iter.max = 1000,nstart = 100)
table(km1$cluster)
Tes_Clusters <- km1$cluster

#### Rename clusters
Tes_Clusters[Tes_Clusters == "1"] <- "BLIS"  ; Tes_Clusters[Tes_Clusters == "2"] <- "LAR"
Tes_Clusters[Tes_Clusters == "3"] <- "IM" ; Tes_Clusters[Tes_Clusters == "4"] <- "MES"

#### Outputs
write.csv(Tes_Clusters, "./results/tables/FUSCCTNBC_old_subtypes.csv")

###################--------------new pipeline-----------------##################
# filtered_fpkm_new_pipe_TT
logexpr_new_TT <- log2(filtered_fpkm_new_pipe_TT+1)
Tes_Expr <- logexpr_new_TT[ProteinCoding, ]
dim(Tes_Expr)

### Expression dataset for clustering
Tes_Expr <- Tes_Expr[Clustering_Genes, ] %>% na.omit()
dim(Tes_Expr)
### 1854 genes(consistent with old pipeline)

### clustering by k-means
set.seed(1234)
km1 <- kmeans(t(Tes_Expr),4,iter.max = 1000,nstart = 100)
table(km1$cluster)
Tes_Clusters <- km1$cluster

#### Rename clusters
Tes_Clusters[Tes_Clusters == "1"] <- "IM"  ; Tes_Clusters[Tes_Clusters == "2"] <- "LAR"
Tes_Clusters[Tes_Clusters == "3"] <- "BLIS" ; Tes_Clusters[Tes_Clusters == "4"] <- "MES"

#### Outputs
write.csv(Tes_Clusters, "./results/tables/FUSCCTNBC_new_subtypes.csv")
rm(list=ls())

# ----------------------Subtyping Column Table---------------------- #
#1. data import
FUSCCTNBC_old_subtypes <- read.csv("./results/tables/FUSCCTNBC_old_subtypes.csv")
colnames(FUSCCTNBC_old_subtypes) <- c("patient","subtype_old")
table(FUSCCTNBC_old_subtypes$subtype_old)
FUSCCTNBC_new_subtypes <- read.csv("./results/tables/FUSCCTNBC_new_subtypes.csv")
colnames(FUSCCTNBC_new_subtypes) <- c("patient","subtype_new")
table(FUSCCTNBC_new_subtypes$subtype_new)
FUSCCTNBC_new_subtypes$patient <- gsub("_rep","",FUSCCTNBC_new_subtypes$patient)
FUSCCTNBC_subtypes <- merge(FUSCCTNBC_old_subtypes,FUSCCTNBC_new_subtypes,by="patient")
rm(FUSCCTNBC_old_subtypes,FUSCCTNBC_new_subtypes)
write.csv(FUSCCTNBC_subtypes,"./results/tables/FUSCCTNBC_subtypes.csv")

# 2. Compare the results
## check the subtype status
result = 0
patients = c()
for(i in 1:length(FUSCCTNBC_subtypes$patient)){
  temp <- sum(FUSCCTNBC_subtypes$subtype_old[i] %in% FUSCCTNBC_subtypes$subtype_new[i])
  if(temp == 1){
    check <- FUSCCTNBC_subtypes$patient[i]
    patients <- append(patients,check)
  }
  result = result+temp
}
### 340 patients matched; 20 patients unmatched

df <- xtabs(~ subtype_old + subtype_new, data = FUSCCTNBC_subtypes)
tabel1 <- as.data.frame(df) %>% melt()
library(plyr)
tabel1 <- ddply(tabel1, "subtype_old", transform,
                percent_old = value / sum(value) * 100)
tabel1 <- ddply(tabel1, "subtype_new", transform,
                percent_new = value / sum(value) * 100)

## plot 
library(ggplot2)
lev <- c("BLIS","IM","MES","LAR")
tabel1$subtype_old <- factor(tabel1$subtype_old,levels = rev(lev))
tabel1$subtype_new <- factor(tabel1$subtype_new,levels = lev)
p1 <- ggplot(tabel1, aes(x=subtype_old, y=percent_old/100, fill=subtype_new)) +
  geom_bar(stat="identity", colour="black")+
  guides(fill=guide_legend(reverse=TRUE))+
  scale_fill_manual(values = c("#ef4922","#7bc242","#3baade","#8f80ba"))+theme_classic()+ 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size=12, face="bold"),
        axis.title.x = element_text(angle = 0, hjust = 0.5, size=14, face="bold"),
        panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        legend.position="none",
        plot.margin=unit(rep(1,4),'lines')) + 
  scale_y_continuous(expand = c(0,0),position = "right")+
  ylab("Proportion")+coord_flip();p1
# +scale_y_reverse()
ggsave("./results/charts/subtype_old_new.png",p1,width = 4,height = 4,dpi = 300,bg = "transparent")


tabel1$subtype_old <- factor(tabel1$subtype_old,levels = rev(lev))
tabel1$subtype_new <- factor(tabel1$subtype_new,levels = rev(lev))
p2 <- ggplot(tabel1, aes(x=subtype_new, y=percent_new/100, fill=subtype_old)) +
  geom_bar(stat="identity", colour="black")+
  scale_fill_manual(values = c("#ef4922","#7bc242","#3baade","#8f80ba"))+theme_classic()+ 
  ylab("Proportion")+
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, size=12, face="bold"),
        axis.title.y = element_text(angle = 90, hjust = 0.5, size=14, face="bold"),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        legend.position="none")+
  scale_y_continuous(expand = c(0,0),position = "left");p2

ggsave("./results/charts/subtype_new_old.png",p2,width = 4,height = 4,dpi = 300,bg = "transparent")

## remove all
rm(list=ls())

