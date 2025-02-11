## R script to compute assay measures from seahorse data
## Script written by Dr. Ankita, Bioinformatician, ACHRI, University of Calgary

### install required
list.of.packages <- c("reshape2", "dplyr", "tidyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#setwd("C:/Users/lab/Desktop/Seahorse_script")

library(dplyr)
library(reshape2)
library(tidyr)
require(purrr)
library(tidyverse)

var <- readline("Enter inputfile : ");
x<-read.table(var, sep="\t", header = TRUE);
start.time <- Sys.time()

Basal_repiration<-x %>% group_by(Group, Measurement) %>% filter( Measurement %in% (3:7))%>% mutate(row = row_number()) %>% tidyr::pivot_wider(names_from = Group, values_from = OCR) %>% select(-row)%>% as.data.frame()%>% mutate(Function="Basal respiration") %>% as.data.frame()
Oligomycin<-x %>% group_by(Group, Measurement) %>% filter( Measurement %in% (12:14))%>% mutate(row = row_number()) %>% tidyr::pivot_wider(names_from = Group, values_from = OCR) %>% select(-row)%>% as.data.frame()%>% mutate(Function="Oligomycin") %>% as.data.frame()
MRC<-x %>% group_by(Group, Measurement) %>% filter( Measurement %in% (16:20))%>% mutate(row = row_number()) %>% tidyr::pivot_wider(names_from = Group, values_from = OCR) %>% select(-row)%>% as.data.frame()%>% mutate(Function="Maximum Respiration Capacity") %>% as.data.frame()
NMR<-x %>% group_by(Group, Measurement) %>% filter( Measurement %in% (25:27))%>% mutate(row = row_number()) %>% tidyr::pivot_wider(names_from = Group, values_from = OCR) %>% select(-row)%>% as.data.frame()%>% mutate(Function="Non-mitochondrial Respiration") %>% as.data.frame()




#target1_ALR <- c(5:7)
#target2_ALR <- c(12:14)


#y<-x %>% group_by(Group) %>% filter (Measurement %in% target) %>% mutate(row = row_number()) %>% tidyr::pivot_wider(names_from = Group, values_from = OCR) %>% select(-row) 


## 
groups<-unique(x$Group)

out_cycle3<-c()
out_cycle4<-c()

out_cycle5<-c()
out_cycle6<-c()
out_cycle7<-c()


out_cycle12<-c()
out_cycle13<-c()
out_cycle14<-c()
out_cycle16<-c()
out_cycle17<-c()
out_cycle18<-c()
out_cycle19<-c()
out_cycle20<-c()


out_cycle25<-c()
out_cycle26<-c()
out_cycle27<-c()

cycle_5_12 <-c()
cycle_6_13 <-c()
cycle_7_14 <-c()

cycle_5_25 <-c()
cycle_6_26 <-c()
cycle_7_27 <-c()

cycle_16_3 <-c()
cycle_17_4 <-c()
cycle_18_5 <-c()
cycle_19_6 <-c()
cycle_20_7 <-c()

cycle_12_25 <-c()
cycle_13_26 <-c()
cycle_14_27 <-c()

#cycle_5_7 <-c()


 for (i in seq_along(groups))
 { 
 print (i)

  
    out_cycle3[[i]]<-x %>% filter (Measurement == "3" & Group %in% groups[[i]])
	  out_cycle4[[i]]<-x %>% filter (Measurement == "4" & Group %in% groups[[i]])
  
  out_cycle5[[i]]<-x %>% filter (Measurement == "5" & Group %in% groups[[i]])
  out_cycle6[[i]]<-x %>% filter (Measurement == "6" & Group %in% groups[[i]])
  out_cycle7[[i]]<-x %>% filter (Measurement == "7" & Group %in% groups[[i]])
  
  
  
  out_cycle12[[i]]<-x %>% filter (Measurement == "12" & Group %in% groups[[i]])
  out_cycle13[[i]]<-x %>% filter (Measurement == "13" & Group %in% groups[[i]])
  out_cycle14[[i]]<-x %>% filter (Measurement == "14" & Group %in% groups[[i]])
  out_cycle16[[i]]<-x %>% filter (Measurement == "16" & Group %in% groups[[i]])
  out_cycle17[[i]]<-x %>% filter (Measurement == "17" & Group %in% groups[[i]])
  out_cycle18[[i]]<-x %>% filter (Measurement == "18" & Group %in% groups[[i]])
  out_cycle19[[i]]<-x %>% filter (Measurement == "19" & Group %in% groups[[i]])
  out_cycle20[[i]]<-x %>% filter (Measurement == "20" & Group %in% groups[[i]])
  out_cycle25[[i]]<-x %>% filter (Measurement == "25" & Group %in% groups[[i]])
  out_cycle26[[i]]<-x %>% filter (Measurement == "26" & Group %in% groups[[i]])
  out_cycle27[[i]]<-x %>% filter (Measurement == "27" & Group %in% groups[[i]])
  
 #print (out_cycle5[[i]])
  #print (out_cycle12[[i]])
  
  #cycle_5_7[[i]] <-x %>% filter (Measurement %in% c("5","6","7") & Group %in% groups[[i]])
  
  cycle_5_12[[i]]<-(out_cycle5[[i]][3]-out_cycle12[[i]][3])
  cycle_6_13[[i]]<-(out_cycle6[[i]][3]-out_cycle13[[i]][3])
  cycle_7_14[[i]]<-(out_cycle7[[i]][3]-out_cycle14[[i]][3])
  
   
  cycle_5_25[[i]]<-(out_cycle5[[i]][3]-out_cycle25[[i]][3])
  cycle_6_26[[i]]<-(out_cycle6[[i]][3]-out_cycle26[[i]][3])
  cycle_7_27[[i]]<-(out_cycle7[[i]][3]-out_cycle27[[i]][3])
  
  cycle_16_3[[i]]<-(out_cycle16[[i]][3]-out_cycle3[[i]][3])
  cycle_17_4[[i]]<-(out_cycle17[[i]][3]-out_cycle4[[i]][3])
  cycle_18_5[[i]]<-(out_cycle18[[i]][3]-out_cycle5[[i]][3])
  cycle_19_6[[i]]<-(out_cycle19[[i]][3]-out_cycle6[[i]][3])
  cycle_20_7[[i]]<-(out_cycle20[[i]][3]-out_cycle7[[i]][3])
  
  
  cycle_12_25[[i]]<-(out_cycle12[[i]][3]-out_cycle25[[i]][3])
  cycle_13_26[[i]]<-(out_cycle13[[i]][3]-out_cycle26[[i]][3])
  cycle_14_27[[i]]<-(out_cycle14[[i]][3]-out_cycle27[[i]][3])
  
  
  
  #print (cycle_5_12[[i]])
  }

## AUtomation: Convert list into dataframe 
y<-NULL;
files <- ls (pattern = "^cycle")
for (i in seq_along(files))
{
df<-get(files[[i]])
L<-list(df)
td<-map_dfr(df, ~unlist(.x) %>% t() %>% as.data.frame())
#td<-map_df(df, ~as.data.frame(t(.)))
new_df<-as.data.frame(t(td))

for (j in seq_along(groups)) { 

colnames(new_df)[j] <- c(groups[j]);
row.names= FALSE;
}

final_df <-new_df  %>% mutate(across(where(is.numeric), ~ round(., 8))) %>% mutate(Measurement=files[i]) %>% mutate(Function = case_when(files[i]=="cycle_5_12" ~ 'ATP-LINKED RESPIRATION', files[i]=="cycle_6_13" ~ 'ATP-LINKED RESPIRATION', files[i]=="cycle_7_14" ~ 'ATP-LINKED RESPIRATION', files[i]=="cycle_5_25" ~ 'TOTAL MITOCHONDRIAL RESPIRATION', files[i]=="cycle_6_26" ~ 'TOTAL MITOCHONDRIAL RESPIRATION', files[i]=="cycle_7_27" ~ 'TOTAL MITOCHONDRIAL RESPIRATION', files[i]=="cycle_16_3" ~ "RESERVE CAPACITY", files[i]=="cycle_17_4" ~ "RESERVE CAPACITY", files[i]=="cycle_18_5" ~ "RESERVE CAPACITY", files[i]=="cycle_19_6" ~ "RESERVE CAPACITY", files[i]=="cycle_20_7" ~ "RESERVE CAPACITY", files[i]=="cycle_12_25" ~ "PROTON LEAKS", files[i]=="cycle_13_26" ~ "PROTON LEAKS", files[i]=="cycle_14_27" ~ "PROTON LEAKS")) 

#df_new<- data.frame(matrix(unlist(df), ncol=length(df))) %>% mutate(across(where(is.numeric), ~ round(., 8))) %>% mutate(Measurement=files[i]) %>% mutate(Function = case_when(files[i]=="cycle_5_12" ~ 'ATP-LINKED RESPIRATION', files[i]=="cycle_6_13" ~ 'ATP-LINKED RESPIRATION', files[i]=="cycle_7_14" ~ 'ATP-LINKED RESPIRATION', files[i]=="cycle_5_25" ~ 'TOTAL MITOCHONDRIAL RESPIRATION', files[i]=="cycle_6_26" ~ 'TOTAL MITOCHONDRIAL RESPIRATION', files[i]=="cycle_7_27" ~ 'TOTAL MITOCHONDRIAL RESPIRATION', files[i]=="cycle_16_3" ~ "RESERVE CAPACITY", files[i]=="cycle_17_4" ~ "RESERVE CAPACITY", files[i]=="cycle_18_5" ~ "RESERVE CAPACITY", files[i]=="cycle_19_6" ~ "RESERVE CAPACITY", files[i]=="cycle_20_7" ~ "RESERVE CAPACITY", files[i]=="cycle_12_25" ~ "PROTON LEAKS", files[i]=="cycle_13_26" ~ "PROTON LEAKS", files[i]=="cycle_14_27" ~ "PROTON LEAKS")) 
#df_new<-t(map_dfr(df, ~as.data.frame(t(.)))) %>% mutate(across(where(is.numeric), ~ round(., 8))) %>% mutate(Measurement=files[i]) %>% mutate(Function = case_when(files[i]=="cycle_5_12" ~ 'ATP-LINKED RESPIRATION', files[i]=="cycle_6_13" ~ 'ATP-LINKED RESPIRATION', files[i]=="cycle_7_14" ~ 'ATP-LINKED RESPIRATION', files[i]=="cycle_5_25" ~ 'TOTAL MITOCHONDRIAL RESPIRATION', files[i]=="cycle_6_26" ~ 'TOTAL MITOCHONDRIAL RESPIRATION', files[i]=="cycle_7_27" ~ 'TOTAL MITOCHONDRIAL RESPIRATION', files[i]=="cycle_16_3" ~ "RESERVE CAPACITY", files[i]=="cycle_17_4" ~ "RESERVE CAPACITY", files[i]=="cycle_18_5" ~ "RESERVE CAPACITY", files[i]=="cycle_19_6" ~ "RESERVE CAPACITY", files[i]=="cycle_20_7" ~ "RESERVE CAPACITY", files[i]=="cycle_12_25" ~ "PROTON LEAKS", files[i]=="cycle_13_26" ~ "PROTON LEAKS", files[i]=="cycle_14_27" ~ "PROTON LEAKS")) 
#assign(paste(files[i],"_df"),df_new)

y<-  rbind(y, final_df)
minus_cycles<-y%>% relocate(Measurement)

}


OUTPUT_SH<-rbind(minus_cycles, Basal_repiration, Oligomycin, NMR, MRC)
write.table(OUTPUT_SH, file="Output_all_measures.txt", sep="\t",quote = FALSE, row.names=FALSE)



end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
print(time.taken)





















