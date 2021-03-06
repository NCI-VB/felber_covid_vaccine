
suppressMessages(library(tidyverse))

setwd("~/Documents/nci_vb_git/felber_covid_vaccine/code")

immune_diff_counts <- read.csv("../results/immune_diff_counts.csv", header = TRUE, check.names = FALSE)
immune_diff_metadata <- read.csv("../results/immune_diff_metadata.csv", header = TRUE, check.names = FALSE)
metadata <- read.csv("../data/metadata.csv", header = TRUE, check.names = FALSE)

row.names(immune_diff_counts) <- immune_diff_counts$Gene
immune_diff_counts$Gene <- NULL

contrast <- c("d23-d22","d2-d1")

test <- contrast[1]
ref <- contrast[2]

samples_in_contrast <- grep(paste(contrast,collapse="|"),colnames(immune_diff_counts), value=TRUE)
df.c <- immune_diff_counts[,samples_in_contrast]

test.samples <- grep(test,samples_in_contrast,value=TRUE)
ref.samples  <- grep(ref, samples_in_contrast,value=TRUE)

df.t <- as.data.frame(t(df.c))
df.t$patient_id <- apply(array(row.names(df.t)),1,function(z) immune_diff_metadata[which(immune_diff_metadata$sample_id == z), "patient_id"])

immune_complex_diff_counts <- as.data.frame(matrix(nrow=length(row.names(df.c)),ncol=length(unique(df.t$patient_id)),dimnames=list(row.names(df.c),paste0(unique(df.t$patient_id),"_",test,"-",ref))))

for(i in 1:nrow(immune_complex_diff_counts)){
  gene <- row.names(immune_complex_diff_counts)[i]
  for(j in 1:ncol(immune_complex_diff_counts)){
    col.nam <- colnames(immune_complex_diff_counts)[j]
    patient_id <- unlist(str_split(col.nam,"_"))[1]
    contrast <- unlist(str_split(col.nam,"_"))[2]
    t.level <- paste(unlist(str_split(contrast,"-"))[1:2],collapse="-")
    r.level <- paste(unlist(str_split(contrast,"-"))[3:4],collapse="-")
    
    test.sample <- paste0(patient_id,"_",t.level)
    ref.sample  <- paste0(patient_id,"_",r.level)
    
    #Do we actually have both levels?
    if(!all(c(test.sample,ref.sample) %in% row.names(df.t))){
      immune_complex_diff_counts[i,j] <- NA
      next
    }
    
    #Are either NA?
    if( is.na(df.t[test.sample,gene]) | is.na(df.t[ref.sample,gene]) ){
      immune_complex_diff_counts[i,j] <- NA
      next                
    }
    
    immune_complex_diff_counts[i,j] <- df.t[test.sample,gene] - df.t[ref.sample,gene]
  }
}

immune_complex_diff_counts <- immune_complex_diff_counts %>% rownames_to_column("Gene")
immune_diff_counts <- immune_diff_counts %>% rownames_to_column("Gene")
merged_diff_immune <- merge(immune_diff_counts,immune_complex_diff_counts,by="Gene")

#Make complex metadata
col.nam <- colnames(immune_complex_diff_counts)[colnames(immune_complex_diff_counts) != "Gene"]
patient_id <- apply(array(col.nam),1,function(z) unlist(str_split(z,"_"))[1]) #  unlist(str_split(col.nam,"_"))[1]
contrast <- rep("vac2_vac1",length(patient_id))
complex_diff_immune_metadata <- data.frame(sample_id=col.nam, patient_id = patient_id, contrast = gsub("-","_",contrast))
combined_metadata <- rbind(immune_diff_metadata,complex_diff_immune_metadata)

merged.df <- metadata %>% select(-sample_id,-timepoint) %>% distinct()
print(head(merged.df))
merged_metadata <- merge(combined_metadata,merged.df, by="patient_id")
merged_metadata <- merged_metadata %>% select(sample_id, everything()) %>% distinct() 

write.csv(merged_diff_immune, "../results/immune_plot_counts.csv", row.names = FALSE, quote = FALSE)
write.csv(merged_metadata, "../results/immune_plot_metadata.csv", row.names = FALSE, quote = FALSE)
