
suppressMessages(library(limma)) #3.38.3
suppressMessages(library(tidyverse)) #1.3.0
suppressMessages(library(edgeR)) #3.24.3
suppressMessages(library(stringr)) #1.4.0

voom_msd <- function (counts, design = NULL, ...) 
{
  out <- list()
  counts <- as.matrix(counts)
  
  n <- nrow(counts)
  if (n < 2L)
    stop("Need at least two genes to fit a mean-variance trend")
  
  if (is.null(design)) {
    stop("Need to provide a contrast variable")
  }
  
  y <- t(log2(t(counts + 0.5)))
  fit <- lmFit(y, design, ...)
  
  out$E <- y
  out$design <- design
  new("EList", out)
}

df <- read.csv("../results/naive_dataset.csv", header = TRUE)
genenames <- df[["Gene"]]
if(any(make.names(colnames(df))!=colnames(df))){
  print("Error: The following counts matrix column names are not valid:\n")
  print(colnames(df)[make.names(colnames(df))!=colnames(df)])
  
  print("Likely causes are columns starting with numbers or other special characters eg spaces.")
  stop("Bad column names.")
}

samples_for_deg_analysis = colnames(df)
samples_for_deg_analysis <- samples_for_deg_analysis[samples_for_deg_analysis != "Gene"]

df.m <- df[,samples_for_deg_analysis]

targetfile <- read.csv("../data/metadata.csv", header=TRUE)
targetfile <- targetfile[match(colnames(df.m),targetfile$sample_id),]
targetfile <- targetfile[rowSums(is.na(targetfile)) != ncol(targetfile), ]
df.m <- df.m[,match(targetfile$sample_id,colnames(df.m))]

x <- 2^df.m

ordered_covariates=c("timepoint","patient_id")

ordered_covariates=ordered_covariates[order(ordered_covariates!="timepoint")]

targetfile <- targetfile %>% select(sample_id,one_of(ordered_covariates)) %>% as.data.frame()

row.names(targetfile) <- targetfile$sample_id
dm.formula <- as.formula(paste("~0 +", paste(ordered_covariates, sep="+", collapse="+")))
design=model.matrix(dm.formula, targetfile)

colnames(design) <- str_replace_all(colnames(design), "timepoint", "")

v <- voom_msd(x,design=design,plot=FALSE)

rownames(v$E) <- genenames
as.data.frame(v$E) %>% rownames_to_column("Gene") -> df.voom
fit <- lmFit(v, design)

contrasts_of_interest <- c( "d2-d1",
                            "d23-d22",
                            "d8-d1",
                            "(d23-d22)-(d2-d1)")

cm <- makeContrasts(contrasts = contrasts_of_interest, levels=design)

fit2 <- contrasts.fit(fit, cm)

fit2 <- eBayes(fit2)


logFC = fit2$coefficients
colnames(logFC)=paste(gsub("-","_",colnames(logFC)),"logFC",sep="_")
tstat = fit2$t
colnames(tstat)=paste(gsub("-","_",colnames(tstat)),"tstat",sep="_")
FC = 2^fit2$coefficients
FC = ifelse(FC<1,-1/FC,FC)
colnames(FC)=paste(gsub("-","_",colnames(FC)),"FC",sep="_")
pvalall=fit2$p.value
colnames(pvalall)=paste(gsub("-","_",colnames(pvalall)),"pval",sep="_")
pvaladjall=apply(pvalall,2,function(x) p.adjust(x,"BH"))
colnames(pvaladjall)=paste(gsub("-","_",colnames(fit2$coefficients)),"adjpval",sep="_")

finalres=as.data.frame(cbind(FC, logFC, tstat, pvalall, pvaladjall))

finalres %>% rownames_to_column("Gene") -> finalres
print(paste0("Total number of genes included: ", nrow(finalres)))

call_me_alias<-colnames(finalres)
colnames(finalres)<-gsub("\\(|\\)","",call_me_alias)
write.csv(finalres,file="../results/naive_deg_analysis.csv",row.names = FALSE, quote = FALSE)
