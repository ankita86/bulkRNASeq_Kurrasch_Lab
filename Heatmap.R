#Customised heatmap for priotised gene set using log normalised log counts of DeSeq2 

library(dplyr)
library("DESeq2")
library("ggplot2")
library("EnhancedVolcano")
library ("ComplexHeatmap")



### select genes of interest and map with normalised counts
genes_selected <-read.table("down_priortised_genes_WT_DM5a.txt", header=T, sep="\t")

### get normalised counts from The DESeqDataSet (dds) created using script deseq.R

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


			  
