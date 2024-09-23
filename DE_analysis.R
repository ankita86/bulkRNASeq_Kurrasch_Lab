## R code -  Differential gene analysis using DESeq2 and gene ontology analysis -

## Some parts of the code were adapted from Thyme et al., bioRxiv preprint, 2024

library(dplyr)
library("DESeq2")
library("ggplot2")
library ("ComplexHeatmap")

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


#### Customised heatmap for priotised gene set using log normalised log counts of DeSeq2 #######
 
### select genes of interest and map with normalised counts
genes_selected <-read.table("down_priortised_genes_WT_DM5a.txt", header=T, sep="\t")

### get normalised counts from The DESeqDataSet (dds) 

countsdata <-log2(counts(dds, normalized=TRUE) + 1)
gene_name <-read.csv("llgeneid_genename.csv")
	 dataframe_counts <-as.data.frame(countsdata)
	 dataframe_counts$LLgeneID<-row.names(dataframe_counts)
	 
dataframe_res <-inner_join(dataframe_counts,gene_name,by="LLgeneID")
	
	
clean <-dataframe_res[dataframe_res$LLgeneSymbol %in% genes_selected$LLgeneSymbol, ]

## merge counts and genes selected to add category column
clean<-inner_join(clean, genes_selected, by="LLgeneSymbol") 
clean <-clean %>% arrange(Category)


	
rownames(clean) <- clean$LLgeneSymbol
clean<-clean[-c(13:14)]
clean_ordered_cols<-select(clean, WT1, WT2, WT3, "5a1_1", "5a1_2", "5a1_3", DM1, DM2, DM3)

mat_gene<-as.matrix(clean_ordered_cols)


 
col_names <- c("WT_1", "WT_2", "WT_3", 
                           expression(~italic("5a"^{-"/-"})*"_1"), expression(~italic("5a"^{-"/-"})*"_2"), expression(~italic("5a"^{-"/-"})*"_3"),
                           expression(~italic("5a"^{-"/-"} * ";" * "5b"^{-"/-"}) * "_1"), expression(~italic("5a"^{-"/-"} * ";" * "5b"^{-"/-"}) * "_2"), expression(~italic("5a"^{-"/-"} * ";" * "5b"^{-"/-"})* "_3"))
 
 
 

   
   

#### Column annotation: customised blocks of annotation ->
grouping_factor <- factor(rep(c("group1", "group2", "group3"), each = 3), 
                          levels = c("group1", "group2", "group3"))


col_annotation <- HeatmapAnnotation(
  foo = anno_block(
    gp = gpar(fill = c("grey", "pink", "red")),  # Colors for each group
    labels = c("WT", expression(~italic("5a"^{-"/-"})),expression(~italic("5a"^{-"/-"} * ";" * "5b"^{-"/-"})) ),  # Labels for each block
    labels_gp = gpar(col = "black", fontsize = 12, fontface = "bold"),  # White color for visibility
    labels_rot = 0,  # No rotation, keeping the labels horizontal
    labels_just = "center"  # Center the labels within the blocks
  )
)
   
   




annotation_df <- genes_selected
rownames(annotation_df) <- annotation_df$LLgeneSymbol
annotation_df <- annotation_df[match(rownames(mat_gene), annotation_df$LLgeneSymbol), ]

## down -> WT_vs_DM5a
row_annotation <- rowAnnotation(
  Category = annotation_df$Category,
  col = list(Category = c("epilepsy" = "skyblue", 
                          "other neurological disorders" = "yellow", 
                         "regulation of calcium ion import" = "green", 
						 "sleep disorders" = "orange"))
)


## up-> WT_vs_DM5a
#row_annotation <- rowAnnotation(
 # Category = annotation_df$Category,
  #col = list(Category = c("apoptosis" = "skyblue", 
   #                       "inflammation and immune response" = "yellow", 
    #                     "metabolism" = "green"))
					
#)


#for up-
#tiff("Heatmap_up_priortised_genes_WT_DM5a.tiff", width = 9, height=4, units = "in",res = 300)
 
 #for down-
 tiff("Heatmap_down_priortised_genes_WT_DM5a.tiff", width = 9, height=6, units = "in",res = 300)




 draw(Heatmap(
     t(scale(t(mat_gene))), 
     name = "Expression", 
     left_annotation = row_annotation,  
     show_row_names = TRUE, 
     cluster_rows = FALSE, 
     cluster_columns = FALSE, 
     row_names_gp = gpar(fontface = "italic"), 
     column_labels = col_names,  # Set empty strings as column labels to match the matrix
     column_split = grouping_factor,  # Split columns into 3 groups
     heatmap_legend_param = list(direction = "horizontal", legend_width = unit(6, "cm")), 
     top_annotation = col_annotation, column_title = NULL
 ),  heatmap_legend_side = "bottom", legend_grouping = "original")


dev.off()	 


			  




