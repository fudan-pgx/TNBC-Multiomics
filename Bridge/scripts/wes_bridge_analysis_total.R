## Title: TNBC Bridge WES Analysis      ##
## Author: Qingwang Chen                ##
## Date: 2022-03-02                     ##  
## Version: V2.0                        ##
##########################################

################-------------------------Preperation---------------------################
# Set work directory and import libraries
rm(list = ls (all = TRUE))
source("./scripts/libraries.R")
source("./scripts/parameter.R")

# 2. Data preprocess
maf_TNscope_new <- annovarToMaf("./data/old_data_new_pipeline_wes/TNscope.ANNOVAR.merged.txt",
                                Center = NULL,refBuild = "hg38",tsbCol = "Tumor_Sample_Barcode",
                                table = "refGene",ens2hugo = TRUE,
                                basename = "./data/old_data_new_pipeline_wes/TNscope.ANNOVAR.merged",
                                sep = "\t",MAFobj = FALSE,sampleAnno = NULL)

maf_TNseq_new <- annovarToMaf("./data/old_data_new_pipeline_wes/TNseq.ANNOVAR.merged.txt",
                                Center = NULL,refBuild = "hg38",tsbCol = "Tumor_Sample_Barcode",
                                table = "refGene",ens2hugo = TRUE,
                                basename = "./data/old_data_new_pipeline_wes/TNseq.ANNOVAR.merged",
                                sep = "\t",MAFobj = FALSE,sampleAnno = NULL)

maf_VarScan_new <- annovarToMaf("./data/old_data_new_pipeline_wes/VarScan.ANNOVAR.merged.txt",
                              Center = NULL,refBuild = "hg38",tsbCol = "Tumor_Sample_Barcode",
                              table = "refGene",ens2hugo = TRUE,
                              basename = "./data/old_data_new_pipeline_wes/VarScan.ANNOVAR.merged",
                              sep = "\t",MAFobj = FALSE,sampleAnno = NULL)

rm(maf_TNscope_new,maf_TNseq_new,maf_VarScan_new)

# 2. Data import and process
final_new <- readRDS("./data/old_data_new_pipeline_wes/TNBC.Extended_WES_Somatic.rds")
final_old <- fread("./data/old_data_old_pipeline_wes/data_mutations_extended_somatic.maf",sep = "\t")
TNscope_new <- fread("./data/old_data_new_pipeline_wes/TNscope.ANNOVAR.merged.maf",sep = "\t")
TNscope_old <- fread("./data/old_data_old_pipeline_wes/299samples.TNscope.PON.cosmic.TN.CloseClusterEvents.PASS.maf",sep = "\t")
TNseq_new <- fread("./data/old_data_new_pipeline_wes/TNseq.ANNOVAR.merged.maf",sep = "\t")
TNseq_old <- fread("./data/old_data_old_pipeline_wes/299samples.TNseq.PON.cosmic.PASS.maf",sep = "\t")
VarScan_new <- fread("./data/old_data_new_pipeline_wes/VarScan.ANNOVAR.merged.maf",sep = "\t")
VarScan_old <- fread("./data/old_data_old_pipeline_wes/299samples.VarScan.TN.SNP_INDEL.Somatic.hc.maf",sep = "\t")
colnames(VarScan_old) <- colnames(TNseq_old)

## Guaranteed Variant_Classification Consistency
table(final_new$Variant_Classification);table(final_old$Variant_Classification)
VC <- unique(final_old$Variant_Classification)
VC_pick <- function(data_df)
{
  data_for_next_analysis <- data_df[data_df$Variant_Classification %in% VC,]
  return(data_for_next_analysis)
}

var_name <- ls(pattern = "new|old")
for(i in var_name){
  n <- paste0(i,"_vc")
  assign(n,VC_pick(get(i)))
  m <- table(get(i)$Variant_Classification)
  print(m)
}

## Guaranteed VAF >=0.08 , sequencing depth >= 8, sequencing read2 >= 2(for final)
sum(final_old_vc$VAF_BamReadcount < 0.08);sum(final_old_vc$depth_BamReadcount < 8)
sum(final_old_vc$alt_dep_BamReadcount < 2)

final_new_vc <- final_new_vc[final_new_vc$VAF >= 0.08 & 
                               final_new_vc$bam_alt_count >= 2 & 
                               final_new_vc$bam_total_count >= 8,]


## 2.3 export and save the data
ana4next <- var_name <- ls(pattern = "_vc")
# dir.create("./results/tables/Variant_Classification")
for(i in ana4next){
  saveRDS(get(i),paste0("./results/tables/Variant_Classification/",i,".rds"))
  print(i)
}


final_new_A <- data.frame(sample_id=final_new_vc$sample_id,
                          Hugo_Symbol=final_new_vc$Hugo_Symbol)
final_old_A <- data.frame(sample_id=final_old_vc$Tumor_Sample_Barcode,
                          Hugo_Symbol=final_old_vc$Hugo_Symbol)

## 2.1 pick the columns for new data
columns_pick<-function(data_new)
{
  data_for_next_analysis <- data.frame(sample_id=data_new$Tumor_Sample_Barcode,
                                       Hugo_Symbol=data_new$Gene.refGene)
  return(data_for_next_analysis)
}

TNscope_new_A <- columns_pick(TNscope_new_vc) %>% as.data.frame()
TNseq_new_A <- columns_pick(TNseq_new_vc) %>% as.data.frame()
VarScan_new_A <- columns_pick(VarScan_new_vc) %>% as.data.frame()

## check patients
patients_id <- unique(TNscope_new_A$sample_id)
setdiff(patients_id,unique(TNseq_new_A$sample_id));setdiff(patients_id,unique(VarScan_new_A$sample_id))
setdiff(unique(TNseq_new_A$sample_id),patients_id);setdiff(unique(VarScan_new_A$sample_id),patients_id)

### Define a high confidence mutation set (mutation detected in over 2/3 callers)
if (T){
  mutation_set1 <- dplyr::intersect(TNscope_new_A,TNseq_new_A)
  mutation_set2 <- dplyr::intersect(TNscope_new_A,VarScan_new_A)
  mutation_set3 <- dplyr::intersect(TNseq_new_A,VarScan_new_A)
  mutation_set4 <- Reduce(intersect, list(TNscope_new_A,TNseq_new_A,VarScan_new_A))
  mutation_set_new <- Reduce(union, list(mutation_set1,mutation_set2,mutation_set3,mutation_set4))
} 

saveRDS(mutation_set_new,"./results/tables/mutation_set_new_A.rds")

## 2.2 sample_id match (for old data)
### define barcode2sample_id function
sampleid_convert<-function(data_old)
{
  data_old$sample_id <- data_old$Tumor_Sample_Barcode
  data_old$series_id <- sprintf("%03d",as.numeric(gsub("CNV|\\.|T","",data_old$sample_id)))
  data_old$sample_id <- paste0("FUSCCTNBC",data_old$series_id)
  data_for_next_analysis <- data.frame(sample_id=data_old$sample_id,
                                       Hugo_Symbol=data_old$Hugo_Symbol)
  return(data_for_next_analysis)
}

TNscope_old_A <- sampleid_convert(TNscope_old_vc) %>% as.data.frame()
TNseq_old_A <- sampleid_convert(TNseq_old_vc) %>% as.data.frame()
VarScan_old_A <- sampleid_convert(VarScan_old_vc) %>% as.data.frame()

## check patients
setdiff(patients_id,unique(TNscope_old_A$sample_id));setdiff(unique(TNscope_old_A$sample_id),patients_id)

## filter the same patients
TNscope_old_A <- TNscope_old_A[TNscope_old_A$sample_id %in% patients_id,]
TNseq_old_A <- TNseq_old_A[TNseq_old_A$sample_id %in% patients_id,]
VarScan_old_A <- VarScan_old_A[VarScan_old_A$sample_id %in% patients_id,]

setdiff(patients_id,unique(VarScan_old_A$sample_id));setdiff(unique(VarScan_old_A$sample_id),patients_id)

### Define a high confidence mutation set (mutation detected in over 2/3 callers)
if (T){
  mutation_set1 <- dplyr::intersect(TNscope_old_A,TNseq_old_A)
  mutation_set2 <- dplyr::intersect(TNscope_old_A,VarScan_old_A)
  mutation_set3 <- dplyr::intersect(TNseq_old_A,VarScan_old_A)
  mutation_set4 <- Reduce(intersect, list(TNscope_old_A,TNseq_old_A,VarScan_old_A))
  mutation_set_old <- Reduce(union, list(mutation_set1,mutation_set2,mutation_set3,mutation_set4))
} 

write.csv(mutation_set_old,"./results/tables/mutation_set_old.csv")
saveRDS(mutation_set_old,"./results/tables/mutation_set_old_A.rds")

### remove irrelevant variants
rm(list=ls()[-grep("_A",ls())])

## 2.3 export and save the data
ana4next <- ls()
for(i in ana4next){
  write.csv(get(i),paste0("./results/tables/",i,".csv"))
  saveRDS(get(i),paste0("./results/tables/",i,".rds"))
  print(i)
}

# remove all
rm(list=ls())

############---------------------Formal Analysis------------------##############
# 1. Library import
source("./scripts/libraries.R")
source("./scripts/parameter.R")

# 2. Data import
dir <- "./results/tables/"
name <- list.files(dir,pattern = "_A.rds$")

for(i in seq_along(name)){
  filename <- gsub(".rds","",name[i])
  assign(filename,readRDS(paste0(dir,name[i])))
}

# 3. Patients Statistic
## variant numbers consistency
var_name <- ls(pattern = "new|old")
for(i in var_name){
  n <- paste0(i,"_var_num")
  m <- table(get(i)$sample_id) %>% as.data.frame()
  colnames(m) <- c("sample_id",gsub("_A","",i))
  assign(n,m)
}

var_name <- ls(pattern = "_var_num")
FUSCCTNBC_var_num <- Reduce(function(x,y) merge(x,y,by="sample_id",all.x=TRUE),
                            list(TNscope_new_A_var_num,
                                 TNscope_old_A_var_num,
                                 TNseq_new_A_var_num,
                                 TNseq_old_A_var_num,
                                 VarScan_new_A_var_num,
                                 VarScan_old_A_var_num,
                                 mutation_set_new_A_var_num,
                                 mutation_set_old_A_var_num,
                                 final_new_A_var_num,
                                 final_old_A_var_num),accumulate =FALSE)

rownames(FUSCCTNBC_var_num) <- FUSCCTNBC_var_num$sample_id
FUSCCTNBC_var_num$sample_id <- NULL
write.csv(FUSCCTNBC_var_num,"./results/tables/FUSCCTNBC_var_num.csv")

## easy_test
plot(FUSCCTNBC_var_num$TNscope_new,FUSCCTNBC_var_num$TNscope_old)
cor.test(FUSCCTNBC_var_num$mutation_set_new,FUSCCTNBC_var_num$mutation_set_old)
cor.test(FUSCCTNBC_var_num$TNscope_new,FUSCCTNBC_var_num$TNscope_old)
cor.test(FUSCCTNBC_var_num$TNseq_new,FUSCCTNBC_var_num$TNseq_old)
cor.test(FUSCCTNBC_var_num$VarScan_new,FUSCCTNBC_var_num$VarScan_old)
cor.test(FUSCCTNBC_var_num$final_new,FUSCCTNBC_var_num$final_old)

#######The Jaccard similarity index is calculated as:#############
#######Jaccard Similarity = (number of observations in both sets) / (number in either set)
#######Or, written in notation form:
#######J(A, B) = |A∩B| / |A∪B|

c

JI_TNscope <- cal_jaccard(TNscope_new_A,TNscope_old_A,"TNscope")
JI_TNseq <- cal_jaccard(TNseq_new_A,TNseq_old_A,"TNseq")
JI_VarScan <- cal_jaccard(VarScan_new_A,VarScan_old_A,"VarScan")
JI_HC_mutation <- cal_jaccard(mutation_set_new_A,mutation_set_old_A,"Filter Mutation")
JI_FN_mutation <- cal_jaccard(final_new_A,final_old_A,"Final Mutation")
FUSCCTNBC_JI <- Reduce(function(x,y) merge(x,y,by="sample_id",all.x=TRUE),
                       list(JI_TNscope,JI_TNseq,JI_VarScan,JI_HC_mutation,JI_FN_mutation),accumulate =FALSE)
rownames(FUSCCTNBC_JI) <- FUSCCTNBC_JI$sample_id
FUSCCTNBC_JI$sample_id <- NULL
write.csv(FUSCCTNBC_JI,"./results/tables/FUSCCTNBC_JI.csv")

## the consistency of variant types
final_new_vc <- readRDS("results/tables/Variant_Classification/final_new_vc.rds")
final_old_vc <- readRDS("results/tables/Variant_Classification/final_old_vc.rds")
### mutated genes JI
final_new_vc_p <- data.frame(sample_id=final_new_vc$sample_id,
                            Hugo_Symbol=final_new_vc$Hugo_Symbol,
                            Variant_Classification=final_new_vc$Variant_Classification)

final_old_vc_p <- data.frame(sample_id=final_old_vc$Tumor_Sample_Barcode,
                             Hugo_Symbol=final_old_vc$Hugo_Symbol,
                             Variant_Classification=final_old_vc$Variant_Classification)

VC <- unique(final_old_vc_p$Variant_Classification)
for(i in VC)
{
  a <- paste0("final_new_vc_p_A_",i)
  b <- final_new_vc_p[final_new_vc_p$Variant_Classification == i,c(1,2)]
  assign(a,b)
  c <- paste0("final_old_vc_p_A_",i)
  d <- final_old_vc_p[final_old_vc_p$Variant_Classification == i,c(1,2)]
  assign(c,d)
  e <-  paste0("JI_VC_",i)
  f <- cal_jaccard(b,d,i)
  assign(e,f)
}

var=ls(pattern = "JI_VC_")
FUSCCTNBC_VC_JI <- Reduce(function(x,y) merge(x,y,by="sample_id",all.x=TRUE),
                       list(JI_VC_Frame_Shift_Del,JI_VC_Frame_Shift_Ins,JI_VC_In_Frame_Del,  
                            JI_VC_In_Frame_Ins,JI_VC_Missense_Mutation,JI_VC_Nonsense_Mutation,   
                            JI_VC_Nonstop_Mutation,JI_VC_RNA,JI_VC_Silent,
                            JI_VC_Splice_Site,JI_VC_Translation_Start_Site),accumulate =FALSE)
rownames(FUSCCTNBC_VC_JI) <- FUSCCTNBC_VC_JI$sample_id
FUSCCTNBC_VC_JI$sample_id <- NULL
write.csv(FUSCCTNBC_VC_JI,"./results/tables/FUSCCTNBC_VC_JI.csv")

### VAF correlation
final_new_vc_v <- data.frame(sample_id=final_new_vc$sample_id,
                             Hugo_Symbol=final_new_vc$Hugo_Symbol,
                             Chr=final_new_vc$Chromosome,
                             Start_Position=final_new_vc$Start_Position,
                             End_Position=final_new_vc$End_Position,
                             Variant_Classification=final_new_vc$Variant_Classification,
                             VAF=final_new_vc$VAF)
final_new_vc_v$Chr <- gsub("chr","",final_new_vc_v$Chr)

final_old_vc_v <- data.frame(sample_id=final_old_vc$Tumor_Sample_Barcode,
                             Hugo_Symbol=final_old_vc$Hugo_Symbol,
                             Chr=final_old_vc$Chromosome,
                             Start_Position=final_old_vc$Start_Position,
                             End_Position=final_old_vc$End_Position,
                             Variant_Classification=final_old_vc$Variant_Classification,
                             VAF=final_old_vc$VAF_BamReadcount)

#### simple test
final_new_vc_t <- unite(final_new_vc_v,"combine",c("sample_id","Chr","Start_Position","End_Position"), sep="-", remove = F)
final_old_vc_t <- unite(final_old_vc_v,"combine",c("sample_id","Chr","Start_Position","End_Position"), sep="-", remove = F)
final_compare <- merge(final_new_vc_t,final_old_vc_t,by="combine")
cor(final_compare$VAF.x,final_compare$VAF.y,method ="pearson")
# 0.9994512 However, the overlapping set of mutations is too small to do the next step in the analysis.


## Known hotspot cancer-related genes from the article
cancer_genes <- c("TP53","PIK3CA","KMT2C","PTEN","RB1","KMT2D","PIK3R1","FAT3","FOXA1",
                  "ATR","NF1","ATK1","ATRX","NOTCH2","APC","NCOR1","HRAS","ERBB2","CBFB")
percentile <- c(74,18,7,6,4,4,4,4,3,3,3,3,3,3,3,3,2,2,1)
cancer_genes_info <- data.frame(cancer_genes = cancer_genes, percentile=percentile)
write.csv(cancer_genes_info,"./results/tables/cancer_genes_info.csv")

## Statistics on the percentage of occurrence of hotspot genes
gene_mutation_frequency_old <- table(mutation_set_old$Hugo_Symbol) %>% as.data.frame()
hot_gene_mutation_frequency_old<- gene_mutation_frequency_old[gene_mutation_frequency_old$Var1 %in% cancer_genes,]
hot_gene_mutation_info_old <- hot_gene_mutation_frequency_old
hot_gene_mutation_info_old$percentile <- hot_gene_mutation_info_old$Freq/length(unique(mutation_set_old$sample_id))

gene_mutation_frequency_new <- table(mutation_set_new$Hugo_Symbol) %>% as.data.frame()
hot_gene_mutation_frequency_new<- gene_mutation_frequency_new[gene_mutation_frequency_new$Var1 %in% cancer_genes,]
hot_gene_mutation_info_new <- hot_gene_mutation_frequency_new
hot_gene_mutation_info_new$percentile <- hot_gene_mutation_info_new$Freq/length(unique(mutation_set_new$sample_id))
compare_genes <- merge(hot_gene_mutation_info_old,hot_gene_mutation_info_new,by="Var1")

# simple test
plot(compare_genes$percentile.x,compare_genes$percentile.y)
cor.test(compare_genes$percentile.x,compare_genes$percentile.y)

write.csv(compare_genes,"./results/tables/compare_genes_info.csv")

## remove all
rm(list=ls())


# 4. Plot
# Library import
source("./scripts/libraries.R")
source("./scripts/parameter.R")
## the consistency of variant numbers
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
## TNscope
p1 <- ggscatter(outliers_remove(FUSCCTNBC_var_num[,c("TNscope_new","TNscope_old")]),
                x = "TNscope_new", y = "TNscope_old",
                add = "reg.line", conf.int = TRUE, 
                add.params = list(linetype=2, color = "lightblue", fill = "lightgray"))+
  stat_cor(method = "spearman",alternative = "two.sided",r.digits = 4,
           p.digits = 4,label.x = 500, label.y = 200)+
  xlim(0,1000)+ylim(0,1000)+
  labs(title="TNscope", x="", y="")+
  theme(plot.title = element_text(hjust = 0.5));p1

## TNseq
p2 <- ggscatter(outliers_remove(FUSCCTNBC_var_num[,c("TNseq_new","TNseq_old")]),
                x = "TNseq_new", y = "TNseq_old",
                add = "reg.line", conf.int = TRUE, 
                add.params = list(linetype=2, color = "lightblue", fill = "lightgray"))+
  stat_cor(method = "spearman",alternative = "two.sided",r.digits = 4,
           p.digits = 4,label.x = 200, label.y = 100)+
  xlim(0,600)+ylim(0,800)+
  labs(title="TNseq", x="", y="")+
  theme(plot.title = element_text(hjust = 0.5));p2

## VarScan
p3 <- ggscatter(outliers_remove(FUSCCTNBC_var_num[,c("VarScan_new","VarScan_old")]),
                x = "VarScan_new", y = "VarScan_old",
                add = "reg.line", conf.int = TRUE, 
                add.params = list(linetype=2, color = "lightblue", fill = "lightgray"))+
  stat_cor(method = "spearman",alternative = "two.sided",r.digits = 4,
           p.digits = 4,label.x = 200, label.y = 50)+
  xlim(0,500)+ylim(0,500)+
  labs(title="VarScan", x="", y="")+
  theme(plot.title = element_text(hjust = 0.5));p3

## filter Mutation
p4 <- ggscatter(outliers_remove(FUSCCTNBC_var_num[,c("mutation_set_new","mutation_set_old")]),
                x = "mutation_set_new", y = "mutation_set_old",
                add = "reg.line", conf.int = TRUE, 
                add.params = list(linetype=2, color = "lightblue", fill = "lightgray"))+
  stat_cor(method = "spearman",alternative = "two.sided",r.digits = 4,
           p.digits = 4,label.x = 200, label.y = 50)+
  xlim(0,500)+ylim(0,600)+
  labs(title="Filter Mutation", x="", y="")+
  theme(plot.title = element_text(hjust = 0.5));p4

## Final Mutation
p5 <- ggscatter(outliers_remove(FUSCCTNBC_var_num[,c("final_new","final_old")]),
                x = "final_new", y = "final_old",
                add = "reg.line", conf.int = TRUE, 
                add.params = list(linetype=2, color = "lightblue", fill = "lightgray"))+
  stat_cor(method = "spearman",alternative = "two.sided",r.digits = 4,
           p.digits = 4,label.x = 200, label.y = 50)+
  xlim(0,500)+ylim(0,500)+
  labs(title="Final Mutation", x="", y="")+
  theme(plot.title = element_text(hjust = 0.5));p5

p6 <- p1+p2+p3+p4+p5;p6

ggsave('./results/charts/FUSCCTNBC_var_num_single_scatter_merged.png',p6,width=12,height=12)

## the consistency of variant genes
FUSCCTNBC_JI <- read.csv("./results/tables/FUSCCTNBC_JI.csv",row.names = 1)
FUSCCTNBC_JI_forplot <- melt(FUSCCTNBC_JI)
FUSCCTNBC_JI_forplot$variable <- factor(FUSCCTNBC_JI_forplot$variable,
                                        levels = c("VarScan","TNscope","TNseq",
                                                   "Filter.Mutation","Final.Mutation"))

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
                          legend.position = "none")+
                    ylim(0,1);FUSCCTNBC_JI_plo
ggsave('./results/charts/FUSCCTNBC_JI_boxplot.png',FUSCCTNBC_JI_plo,width=8,height=6)

FUSCCTNBC_JI_lineplot <- ggline(FUSCCTNBC_JI_forplot, 
                                x = "variable", 
                                y = "value", 
                                color = "variable",
                                add = c("mean_sd", "jitter","violin"),
                                palette = "jco")+
  labs(title="", x="", y="Jaccard Index")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")+
  ylim(0,1);FUSCCTNBC_JI_lineplot

ggsave('./results/charts/FUSCCTNBC_JI_lineplot.png',FUSCCTNBC_JI_lineplot,width=8,height=6)

## the consistency of variant types
JI_VC_mutation <- read.csv("./results/tables/FUSCCTNBC_VC_JI.csv",row.names = 1)
JI_VC_mutation_forplot <- melt(JI_VC_mutation)

FUSCCTNBC_JI_VC_plo <- ggplot(JI_VC_mutation_forplot,aes(x=variable, y=value,
                                                    color=variable,fill=variable))+
  geom_jitter(alpha = 0.5,size = 0.8)+
  geom_boxplot(alpha = 0.5,width = 0.8,
               outlier.shape = NA)+
  theme_bw()+mytheme+
  scale_colour_manual(values = pal_20[-(1:5)]) +
  scale_fill_manual(values = pal_20[-(1:5)]) +
  labs(title="", x="", y="Jaccard Index")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  ylim(0,1);FUSCCTNBC_JI_VC_plo
ggsave('./results/charts/FUSCCTNBC_JI_VC_boxplot.png',FUSCCTNBC_JI_VC_plo,width=8,height=6)

## the consistency of hotspot mutated genes
compare_genes_info <- read.csv("./results/tables/compare_genes_info.csv",row.names = 1)

p7 <- ggscatter(compare_genes_info, x = "percentile.x", y = "percentile.y",
                add = "reg.line", conf.int = TRUE, 
                add.params = list(linetype=2, color = "lightblue", fill = "lightgray"))+
  stat_cor(method = "pearson",alternative = "two.sided",r.digits = 4,
           p.digits = 4,label.x = 0.4, label.y = 0.2)+
  labs(title="", x="Old pipeline", y="New pipeline")+
  theme(plot.title = element_text(hjust = 0.5))+mytheme+
  geom_label_repel(aes(percentile.x,percentile.y,label=Var1));p7
ggsave('./results/charts/hotspot_genes_precentile_scatter.png',p7,width=4,height=4)

## remove all
rm(list=ls())
