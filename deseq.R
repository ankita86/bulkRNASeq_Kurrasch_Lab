## R code -  Differential gene analysis using DESeq2 and gene ontology analysis -

## Some codes are adapted from Thyme et al., ---

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



## Different contrasts - 
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


DM_WT<-DE(R_DM_WT)
SM_5a_WT<-DE(R_SM_5a_WT)
SM_5b_WT<-DE(R_SM_5b_WT)
SM_5a_DM<-DE(R_SM_5a_DM)
SM_5b_DM<-DE(R_SM_5b_DM)







## Enhanced volcano

comp1_genes <- read.table("complete_DE_files.txt", sep="\t", header=T) # Input is full DE results file for selected contrast.
genes_selected=c(DM_WT$up$LLgeneSymbol[1:10], DM_WT$down$LLgeneSymbol[1:10])

## Top 10 genes - up and down both
EnhancedVolcano(comp1_genes, x = 'log2FoldChange', lab =comp1_genes$LLgeneSymbol, y='padj',title= "DM vs WT", legendPosition = 'right', selectLab=genes_selected, drawConnectors = TRUE)






### Gene Ontology Analysis #########

library("BiocFileCache")
library("org.Dr.eg.db")
library("dbplyr")
library(biomaRt)
library(clusterProfiler)




{
  listMarts()
  ensembl=useMart("ensembl")
  listDatasets(ensembl)
  ensembl=useDataset("drerio_gene_ensembl", mart=ensembl)
  listAttributes(mart=ensembl)
  
  annoDRerio <- getBM(attributes = c("ensembl_gene_id", "ensembl_gene_id_version","entrezgene_id", 
                                     "zfin_id_symbol", "description", "external_gene_name"), mart=ensembl)
}

## function for enrichment analysis
enrich <- function(var_DE)
{

name <- rlang::ensym(var_DE)
name <-as.character(name)
 
 print (name)
 
## filter genes based on log value
up<-var_DE$p_val_0.05[var_DE$p_val_0.05$log2FoldChange>= 1.5,]
down<-var_DE$p_val_0.05[var_DE$p_val_0.05$log2FoldChange <= -1.5,]

## map data with annoDRerio database
anno_gene_list_down <- inner_join(down,annoDRerio,by=c("LLgeneSymbol" ="external_gene_name"))
anno_gene_list_up <- inner_join(up,annoDRerio,by=c("LLgeneSymbol" ="external_gene_name"))

## enrichKEGG function and setReadable to convert gene IDs to gene names
keggPA_down = enrichKEGG(anno_gene_list_down$entrezgene_id, organism="dre",pvalueCutoff=0.05)
kegg_down <- as.data.frame(setReadable(keggPA_down, OrgDb = org.Dr.eg.db, keyType="ENTREZID"))


keggPA_up = enrichKEGG(gene=anno_gene_list_up$entrezgene_id, organism="dre",pvalueCutoff=0.05)
kegg_up <- as.data.frame(setReadable(keggPA_up, OrgDb = org.Dr.eg.db, keyType="ENTREZID"))

## GO enrichment
ego_down <- enrichGO(gene          = anno_gene_list_down$entrezgene_id,
                universe      = names(anno_gene_list_down$LLgeneSymbol),
                OrgDb         = org.Dr.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

ego_up <- enrichGO(gene          = anno_gene_list_up$entrezgene_id,
                universe      = names(anno_gene_list_up$LLgeneSymbol),
                OrgDb         = org.Dr.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)	
write.csv(kegg_down,paste0(name,"_Pathway_down.csv"))
write.csv(kegg_up,paste0(name,"_Pathway_up.csv"))
write.csv(ego_down,paste0(name,"_BP_ontology_down.csv"))
write.csv(ego_up,paste0(name,"_BP_ontology_up.csv"))

#write.csv(kegg_down,file=paste0("Pathway","_down_pathways.csv"))
#write.csv(kegg_up,file=paste0("Pathway","_up_pathways.csv"))
#write.csv(ego_down,file=paste0("BP_ontology","_down_pathways.csv"))
#write.csv(ego_up,file=paste0("BP_ontlogy","_up_pathways.csv"))


##pathway graph

kegg_down$NegLogPAdj <-log10(kegg_down$p.adjust)

kegg_up$NegLogPAdj <--log10(kegg_up$p.adjust)

kegg_forplot <-rbind(kegg_down,kegg_up) ##only use this for making a plot because I flipped the log adjusted p values in order to force the down regulated pathways to be on the left....

return(
          list(Pathway_down=kegg_down,
          Pathway_up=kegg_up, Pathway_all=
kegg_forplot, GO_up=ego_up, GO_down=ego_down ))
	

}


## call function
comp1<-enrich(DM_WT)
comp2 <- enrich(SM_5a_WT)
comp3 <- enrich(SM_5b_WT)


## after calling function - create plots

color<-ifelse(comp1$Pathway_all$NegLogPAdj<0,"blue","yellow")
jpeg(file=paste0("comp1","_pathways.jpeg"))

ggplot(comp1$Pathway_all, aes(y=reorder(Description,NegLogPAdj), x=NegLogPAdj))+geom_bar(stat="identity",fill=color)+
  labs(x="log 10 padj",y="pathway")
  
  
dev.off()


jpeg(file=paste0("comp1","_pathways_up.jpeg"))
ggplot(comp1$Pathway_up, aes(y=Description, x=NegLogPAdj))+geom_bar(stat="identity")
dev.off()

jpeg(file=paste0("comp1","_pathways_down.jpeg"))
ggplot(comp1$Pathway_down, aes(y=Description, x=-NegLogPAdj))+geom_bar(stat="identity")
dev.off()

jpeg(file=paste0("comp1","_GO_up.jpeg"))
dotplot(comp1$GO_up)	
dev.off()

jpeg(file=paste0("comp1","_GO_down.jpeg"))
dotplot(comp1$GO_down)	
dev.off()


