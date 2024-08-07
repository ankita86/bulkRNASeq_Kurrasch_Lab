
 ### Heatmap for priortised genes -

library(dplyr)
library("DESeq2")
library("ggplot2")
library("EnhancedVolcano")
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


### select genes of interest and map with normalised counts
genes_selected <-read.table("list_priortised_genes", header=T, sep="\t")

### Extract gene counts by mapping gene names
countsdata <-log2(counts(dds, normalized=TRUE) + 1)
gene_name <-read.csv("llgeneid_genename.csv")
dataframe_counts <-as.data.frame(countsdata)
dataframe_counts$LLgeneID<-row.names(dataframe_counts)
	 
dataframe_res <-inner_join(dataframe_counts,gene_name,by="LLgeneID")
	
clean <-dataframe_res[dataframe_res$LLgeneSymbol %in% genes_selected$LLgeneSymbol, ]

## merge counts and genes selected to add category column
clean<-inner_join(clean, genes_selected, by="LLgeneSymbol") 
clean <-clean %>% arrange(Category)

## select and order columns of interest to plot -
rownames(clean) <- clean$LLgeneSymbol
clean<-clean[-c(13:14)]
clean_ordered_cols<-select(clean, WT1, WT2, WT3, "5a1_1", "5a1_2", "5a1_3", DM1, DM2, DM3)


## convert dataframe to matrix
mat_gene<-as.matrix(clean_ordered_cols)

## Customising gene column names -
col_names <- c("WT_1", "WT_2", "WT_3", 
                           expression(~italic("5a"^{-"/-"})*"_1"), expression(~italic("5a"^{-"/-"})*"_2"), expression(~italic("5a"^{-"/-"})*"_3"),
                           expression(~italic("5a"^{-"/-"} * ";" * "5b"^{-"/-"}) * "_1"), expression(~italic("5a"^{-"/-"} * ";" * "5b"^{-"/-"}) * "_2"), expression(~italic("5a"^{-"/-"} * ";" * "5b"^{-"/-"})* "_3"))
 
 
 

## Customising column annotation and colour for ComplexHeatmap
 
 
# Define the conditions
condition_colors <- c('WT' = 'grey', 
                      '5a' = 'pink', 
                      'DM' = 'red')

annotation <- data.frame(
  samples = colnames(mat_gene),
  Condition = c('WT', 'WT', 'WT', 
                '5a', '5a', '5a', 
                'DM', 'DM', 'DM')
)

# Create the top annotation
col_annotation <- HeatmapAnnotation(
    Condition = annotation$Condition,
    col = list(Condition = condition_colors), annotation_legend_param = list(
    Condition = list(
      labels = c("WT", expression(~italic("5a"^{-"/-"})),expression(~italic("5a"^{-"/-"} * ";" * "5b"^{-"/-"})) ),
      at = c("WT","5a" , "DM")))
)


#colours <- list('Condition' = c('WT' = 'red2', 'expression("5a"^{-"/-"}' = 'royalBlue', 'DM'="grey"))
 
#col_ann<-HeatmapAnnotation(df = ann, col=colours, which="col")




annotation_df <- genes_selected
rownames(annotation_df) <- annotation_df$LLgeneSymbol
annotation_df <- annotation_df[match(rownames(mat_gene), annotation_df$LLgeneSymbol), ]


## Row annotations -
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


To save high resolution tiff image -
tiff("Heatmap_down_priortised_genes_WT_DM5a.tiff", width = 8, height=4, units = "in",res = 300)
 
 #for down-
 tiff("Heatmap_down_priortised_genes_WT_DM5a.tiff", width = 8, height=6, units = "in",res = 300)


ht_list<-Heatmap(
    t(scale(t(mat_gene))), 
     name = "Expression", 
     left_annotation = row_annotation,
     show_row_names = TRUE, cluster_rows = FALSE, cluster_columns = FALSE, row_names_gp = gpar(fontface = "italic"), column_labels = col_names, top_annotation=col_annotation, heatmap_legend_param = list(direction = "horizontal", legend_width = unit(6, "cm")))

 draw(ht_list, heatmap_legend_side="bottom", annotation_legend_side="right",
      legend_grouping = "original")
	  
dev.off()	  

