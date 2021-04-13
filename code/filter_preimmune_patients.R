
suppressMessages(library(tidyverse))

df <- read.csv("../data/msd_dataset.csv", header = TRUE)

preimmune <- read.csv("../data/preimmune_patients.csv", header = TRUE)

filter <- grep(paste(preimmune$patient_id,collapse="|"),colnames(df),value=TRUE)
df <- df%>% select(one_of(c("Gene",filter)))

write.csv(df, file="../results/immune_dataset.csv",row.names = FALSE, quote = FALSE)