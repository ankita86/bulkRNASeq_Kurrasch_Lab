#To calculate fluorescence variation measurements (ΔF/F) in each frame of all the neurons, followed by calculating the frequency and amplitude of calcium events (with a threshold value of ΔF/F>2) for each larvae recorded.

## Script written by Dr. Ankita, Bioinformatician, ACHRI, University of Calgary

## Enter input file -
var <- readline("Enter inputfile : ");
x<-read.table(var, sep="\t",  row.names=1, header = TRUE);

fish_name<-gsub(".txt","", var)


#x <-read.table ("Data.txt", header=T, sep="\t", row.names=1)

### calculate value based on (F-Fo)/Fo
criteria <- function(x) { 
     ( (x-min(x))/min(x))
	 }
 
min_col<-apply(x, 2, criteria)
#threshold <- readline("Enter inputfile : ");
write.table(min_col, file=paste0(fish_name, "_minimum_criteria1.txt"), sep="\t", quote = FALSE, col.names = NA)

#min_col <- x %>% mutate_all( funs((.-min(.)) / min(.))) 
#min_col <- x %>% mutate_all( funs((.-min(.)) / min(.))) 
 
### values greater than input threshold
threshold <- readline("Enter Threshold : ");
threshold <- as.numeric(threshold)


threshold_func <- function(x) {  
ifelse(x>threshold, x, NA)
#return(test)
}
condition<-apply(min_col, 2, threshold_func)
write.table(condition, file=paste0(fish_name, "_threshold2_criteria.txt"), sep="\t", quote = FALSE, col.names = NA)

## remove columns that contains only NA values
new_condition<-condition[ , colSums(is.na(condition)) < nrow(condition), drop=FALSE]

### count of non-NA values in columns with both non-NA and NA values
count_per_neuron<-colSums(!is.na(new_condition))
mean_counts_per_fish <- mean(count_per_neuron)

## sum of non-NA values in columns with both non-NA and NA values
sum_per_neuron<-apply(new_condition, 2, sum, na.rm=TRUE)
## average of non-NA values in columns with both non-NA and NA values
average_per_neuron<-apply(new_condition, 2, mean, na.rm=TRUE)

#### final average - sum of averages (calculated in the previous step) / no. of columns with non-NA values
final_average_neurons<-sum(average_per_neuron)/length(count_per_neuron)

table<-rbind(count_per_neuron,mean_counts_per_fish, sum_per_neuron, average_per_neuron, final_average_neurons, new_condition)


write.table(table, file=paste0(fish_name, "_final3.txt"), sep="\t", quote = FALSE, col.names = NA)
