suppressMessages(library(tidyverse)) #1.3.0
suppressMessages(library(stringr)) #1.4.0

setwd("/Users/angelmg/Documents/nci_vb_git/felber_covid_vaccine/code")

df.orig <- read.csv("../results/naive_diff_counts.csv", header = TRUE, check.names = FALSE)
annot <- read.csv("../results/naive_diff_metadata.csv", header = TRUE)

row.names(df.orig) <- df.orig$Gene
df.orig$Gene <- NULL

contrast <- c("d23-d22","d2-d1")

test <- contrast[1]
ref <- contrast[2]

samples_in_contrast <- grep(paste(contrast,collapse="|"),colnames(df.orig), value=TRUE)
df.c <- df.orig[,samples_in_contrast]

test.samples <- grep(test,samples_in_contrast,value=TRUE)
ref.samples  <- grep(ref, samples_in_contrast,value=TRUE)

df.t <- as.data.frame(t(df.c))
df.t$patient_id <- apply(array(row.names(df.t)),1,function(z) annot[which(annot$sample_id == z), "patient_id"])

df.ret <- as.data.frame(matrix(nrow=length(row.names(df.c)),ncol=length(unique(df.t$patient_id)),dimnames=list(row.names(df.c),paste0(unique(df.t$patient_id),"_",test,"-",ref))))

for(i in 1:nrow(df.ret)){
  gene <- row.names(df.ret)[i]
  for(j in 1:ncol(df.ret)){
    col.nam <- colnames(df.ret)[j]
    patient_id <- unlist(str_split(col.nam,"_"))[1]
    contrast <- unlist(str_split(col.nam,"_"))[2]
    t.level <- paste(unlist(str_split(contrast,"-"))[1:2],collapse="-")
    r.level <- paste(unlist(str_split(contrast,"-"))[3:4],collapse="-")
    
    test.sample <- paste0(patient_id,"_",t.level)
    ref.sample  <- paste0(patient_id,"_",r.level)
    
    #Do we actually have both levels?
    if(!all(c(test.sample,ref.sample) %in% row.names(df.t))){
      df.ret[i,j] <- NA
      next
    }
    
    #Are either NA?
    if( is.na(df.t[test.sample,gene]) | is.na(df.t[ref.sample,gene]) ){
      df.ret[i,j] <- NA
      next                
    }
    
    df.ret[i,j] <- df.t[test.sample,gene] - df.t[ref.sample,gene]
  }
}

#make diff diff annotation
diff_metadata <- read.csv("../results/naive_diff_metadata.csv", header = TRUE)
col.nam <- colnames(df.ret)[colnames(df.ret) != "Gene"]
patient_id <- apply(array(col.nam),1,function(z) unlist(str_split(z,"_"))[1]) #  unlist(str_split(col.nam,"_"))[1]
contrast <- rep("vac2_vac1",length(patient_id)) #apply(array(col.nam),1,function(z) unlist(str_split(z,"_"))[2])
#t.level <- apply(array(contrast),1,function(z) unlist(str_split(z,"_"))[1]) #unlist(str_split(contrast,"-"))[1]
#r.level <- apply(array(contrast),1,function(z) unlist(str_split(z,"_"))[2]) #unlist(str_split(contrast,"-"))[2]

diff_diff_metadata <- data.frame(sample_id=col.nam, patient_id = patient_id, contrast = gsub("-","_",contrast))
complex_metadata <- rbind(diff_metadata,diff_diff_metadata)

metadata <- read.csv("../data/metadata.csv", header = TRUE)
metadata <- metadata %>% select(-sample_id,-timepoint) %>% distinct()
naive_plot_metadata <- merge(complex_metadata,metadata, by="patient_id")

naive_plot_counts <- merge(df.orig,df.ret,by="row.names")
colnames(naive_plot_counts)[1] <- "Gene"

write.csv(df.all, "../results/naive_plot_counts.csv", row.names = FALSE, quote = FALSE)
write.csv(naive_plot_metadata, "../results/naive_plot_metadata.csv", row.names = FALSE, quote = FALSE)

