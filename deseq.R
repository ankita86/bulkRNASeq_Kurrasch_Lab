## R code -  Differential gene analysis using DESeq2 and gene ontology analysis -

## Some parts of the code were adapted from Thyme et al., bioRxiv preprint, 2024

library(dplyr)
library("DESeq2")
library("ggplot2")
library("EnhancedVolcano")
library(pheatmap)

sampleTable <- read.csv("samples.csv")

treatments=c("WT", "DM", "SM_5a", "SM_5b")

ddsHTSeq <-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory="counts/",
                                       design=~condition)
									   
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,levels=treatments)

dds <-DESeq(ddsHTSeq)

# filtering might need to change depending on sample number
keep <- rowSums(counts(dds) == 0) < 4
dds <- dds[keep,]

### PCA

rlog<-rlog(dds)
pcaplot<-plotPCA(rlog,"condition")

tiff("PCA.tiff", width = 8, height=6, units = "in",res = 300)
## customising PCA plot
pcaplot + scale_color_manual(name="Group", labels = c("WT", expression(~italic("5a"^{-"/-"} * ";" * "5b"^{-"/-"})), expression(~italic("5a"^{-"/-"})), expression(~italic("5b"^{-"/-"}))), values=c("darkgrey", "red", "pink", "skyblue")) + theme_bw()

dev.off()



## Different contrasts ( DM - Double mutant, SM_5a and SM_5b - Single mutatnt 5a and 5b, respectively and WT - Wild type) - 
R_DM_WT <-results(dds,contrast=c("condition","DM","WT"))
R_SM_5a_WT <-results(dds,contrast=c("condition","SM_5a","WT"))
R_SM_5b_WT <-results(dds,contrast=c("condition","SM_5b","WT"))
R_SM_5a_DM <-results(dds,contrast=c("condition","SM_5a","DM"))
R_SM_5b_DM <-results(dds,contrast=c("condition","SM_5b","DM"))

R_DM_WT<-R_DM_WT[order(R_DM_WT$pvalue),]
R_SM_5a_WT<-R_SM_5a_WT[order(R_SM_5a_WT$pvalue),]
R_SM_5b_WT<-R_SM_5b_WT[order(R_SM_5b_WT$pvalue),]
R_SM_5a_DM<-R_SM_5a_DM[order(R_SM_5a_DM$pvalue),]
R_SM_5b_DM<-R_SM_5b_DM[order(R_SM_5b_DM$pvalue),]

## DE function to parse DE analysis results as a dataframe and mapping gene names -

DE <- function(res)
{
    gene_name <-read.csv("llgeneid_genename.csv")
    dataframe_res <-as.data.frame(res)
    dataframe_res$LLgeneID<-row.names(dataframe_res)
    dataframe_res<-dataframe_res[c(7,1:6)]
    res_gene <-inner_join(dataframe_res,gene_name,by="LLgeneID")
    res_gene<-res_gene[!is.na(res_gene$padj),]
    res_gene_p<-subset(res_gene,padj<0.05)   
    res_gene_p_up<-subset(res_gene_p,log2FoldChange>0)
    res_gene_p_down<-subset(res_gene_p,log2FoldChange<0)
    return(list(all=res_gene,p_val_0.05=res_gene_p, up=res_gene_p_up, down=res_gene_p_down))
}

### calling functions to create a list of up-regulated and down-regulated genes with p-value <0.05
DM_WT<-DE(R_DM_WT)
SM_5a_WT<-DE(R_SM_5a_WT)
SM_5b_WT<-DE(R_SM_5b_WT)
SM_5a_DM<-DE(R_SM_5a_DM)
SM_5b_DM<-DE(R_SM_5b_DM)




