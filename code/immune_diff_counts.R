
suppressMessages(library("tidyverse"))

setwd("~/Documents/nci_vb_git/felber_covid_vaccine/code")

df <- read.csv("../results/immune_dataset.csv", header = TRUE,  check.names = FALSE)
annot <- read.csv("../data/metadata.csv", header = TRUE)

annot <- annot %>% filter(sample_id %in% colnames(df))

var_of_interest <- "timepoint"
contrasts <- c("d2-d1","d8-d1","d22-d1","d23-d22")

for(n in seq_along(contrasts)){
  c <- contrasts[n]
  print(c)
  test <- unlist(str_split(c,"-"))[1]
  ref  <- unlist(str_split(c,"-"))[2]
  
  annot.c <- annot %>% filter(!!rlang::sym(var_of_interest) %in% c(test,ref))
  samples_in_contrast <- annot.c$sample_id
  
  df.c <- df[ , c("Gene",samples_in_contrast)]
  row.names(df.c) <- df$Gene
  df.c$Gene <- NULL
  
  test.samples <- annot.c %>% filter(!!rlang::sym(var_of_interest) %in% test) %>% pull(sample_id)
  ref.samples  <- annot.c %>% filter(!!rlang::sym(var_of_interest) %in% ref)  %>% pull(sample_id)
  
  df.t <- as.data.frame(t(df.c))
  
  df.t$patient_id <- apply(array(row.names(df.t)),1,function(z) annot[which(annot$sample_id == z), "patient_id"])
  
  df.ret <- as.data.frame(matrix(nrow=length(row.names(df.c)),ncol=length(unique(df.t$patient_id)),dimnames=list(row.names(df.c),paste0(unique(df.t$patient_id),"_",test,"-",ref))))
  
  for(i in 1:nrow(df.ret)){
    gene <- row.names(df.ret)[i]
    for(j in 1:ncol(df.ret)){
      col.nam <- colnames(df.ret)[j]
      patient_id <- unlist(str_split(col.nam,"_"))[1]
      contrast <- unlist(str_split(col.nam,"_"))[2]
      t.level <- unlist(str_split(contrast,"-"))[1]
      r.level <- unlist(str_split(contrast,"-"))[2]
      
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
  
  df.ret <- df.ret %>% rownames_to_column("Gene")
  if(n == 1){
    immune_diff_counts <- df.ret
  }else{
    immune_diff_counts <- merge(immune_diff_counts,df.ret,by="Gene")
  }
}

col.nam <- colnames(immune_diff_counts)[colnames(immune_diff_counts) != "Gene"]
patient_id <- apply(array(col.nam),1,function(z) unlist(str_split(z,"_"))[1])
contrast <- apply(array(col.nam),1,function(z) unlist(str_split(z,"_"))[2])
immune_diff_metadata <- data.frame(sample_id=col.nam, patient_id = patient_id, contrast = gsub("-","_",contrast))



write.csv(immune_diff_counts, "../results/immune_diff_counts.csv", row.names = FALSE, quote = FALSE)
write.csv(immune_diff_metadata, "../results/immune_diff_metadata.csv", row.names = FALSE, quote = FALSE)
