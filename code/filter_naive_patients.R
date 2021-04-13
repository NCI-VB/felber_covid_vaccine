
suppressMessages(library(tidyverse)) #v1.3.0

setwd("/Users/angelmg/Documents/nci_vb_git/felber_covid_vaccine/code")

df <-  read.csv("../data/msd_dataset.csv", header=TRUE)
preimmune <- read.csv("../data/preimmune_patients.csv", header=TRUE)

filter <- grep(paste(preimmune$patient_id,collapse="|"),colnames(df),value=TRUE)
df <- df%>% select(-one_of(filter))

write.csv(df,file="../results/naive_dataset.csv",row.names = FALSE, quote = FALSE)
