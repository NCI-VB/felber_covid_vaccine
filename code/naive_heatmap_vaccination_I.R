require(pheatmap)
require(dendsort)

suppressMessages(library(colorspace))
suppressMessages(library(dendsort))
suppressMessages(library(pheatmap))
suppressMessages(library(tidyverse))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dplyr))

doheatmap <- function(dat, clus, clus2, ht, rn, cn, col,scale_range) {
  
  col_filter <- apply(as.data.frame(dat), 2, function(z) all(is.na(z)))
  row_filter <- apply(as.data.frame(dat), 1, function(z) length(z[is.na(z)]) > 10)
  dat <- dat[!row_filter, !col_filter]
  
  col.pal <- diverging_hcl(n=100, palette=col)
  if (FALSE) {
    col.pal = rev(col.pal)
  }
  # define metrics for clustering
  drows1 <- "correlation"
  dcols1 <- "correlation"
  
  minx = floor(min(dat))
  maxx = ceiling(max(dat))
  
  
  breaks = seq(-1*scale_range, scale_range, length=100)
  legbreaks = seq(-1*scale_range, scale_range, length=5)
    
  breaks = sapply(breaks, signif, 4)
  legbreaks = sapply(legbreaks, signif, 4)
  
  treeheight <- 25
  
  hm.parameters <- list(
    dat, 
    color=col.pal,
    legend_breaks=legbreaks,
    cellwidth=if(0 > 0){ 14 }else{ NA }, 
    cellheight=if(0 > 0){ 14 }else{ NA }, 
    scale="none",
    treeheight_col=treeheight,
    treeheight_row=treeheight,
    kmeans_k=NA,
    breaks=breaks,
    # height=80,
    fontsize=14,
    fontsize_row=if(0 > 0){ NA }else{ 14 },
    fontsize_col=if(0 > 0){ NA }else{ 14 },
    show_rownames=rn, 
    show_colnames=cn,
    cluster_rows=FALSE, 
    cluster_cols=FALSE,
    cutree_rows=1,
    annotation_col = annotation_col,
    annotation_colors = annot_col,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    labels_col = labels_col,
    annotation_names_col = FALSE,
    na_col = "#000000"
  )
  
  # print('calculated mat')
  
  p <- do.call("pheatmap", c(hm.parameters))
  return(p)
}

df <- read.csv("../results/naive_plot_counts.csv", header = TRUE, check.names = FALSE)

samples_to_include = c("V002_d2-d1","V002_d22-d1","V002_d8-d1","V004_d22-d1","V004_d8-d1","V004_d2-d1","V005_d22-d1","V005_d8-d1","V005_d2-d1","V006_d22-d1","V006_d8-d1","V006_d2-d1","V007_d8-d1","V007_d2-d1","V007_d22-d1","V008_d8-d1","V008_d2-d1","V008_d22-d1","V009_d22-d1","V009_d8-d1","V009_d2-d1","V010_d22-d1","V010_d8-d1","V010_d2-d1","V011_d8-d1","V011_d22-d1","V011_d2-d1","V013_d22-d1","V013_d2-d1","V013_d8-d1","V015_d22-d1","V015_d2-d1","V015_d8-d1","V016_d22-d1","V016_d8-d1","V016_d2-d1","V017_d8-d1","V017_d22-d1","V017_d2-d1","V019_d22-d1","V019_d2-d1","V019_d8-d1","V020_d2-d1","V020_d8-d1","V020_d22-d1","V021_d8-d1","V021_d22-d1","V021_d2-d1","V022_d8-d1","V022_d22-d1","V022_d2-d1","V024_d8-d1","V024_d22-d1","V024_d2-d1","V025_d8-d1","V025_d22-d1","V025_d2-d1","V026_d8-d1","V026_d22-d1","V026_d2-d1","V027_d8-d1","V027_d22-d1","V027_d2-d1","V028_d8-d1","V028_d22-d1","V028_d2-d1","V029_d2-d1","V029_d8-d1","V029_d22-d1","V030_d8-d1","V030_d22-d1","V030_d2-d1","V031_d22-d1","V031_d8-d1","V031_d2-d1","V033_d22-d1","V033_d8-d1","V033_d2-d1","V034_d8-d1","V034_d22-d1","V034_d2-d1","V036_d2-d1","V036_d8-d1","V036_d22-d1","V037_d22-d1","V037_d8-d1","V037_d2-d1","V039_d22-d1","V039_d8-d1","V039_d2-d1","V041_d22-d1","V041_d2-d1","V041_d8-d1","V043_d2-d1","V043_d8-d1","V043_d22-d1","V045_d22-d1","V045_d2-d1","V045_d8-d1","V048_d22-d1","V048_d2-d1","V048_d8-d1","V049_d22-d1","V049_d2-d1","V049_d8-d1","V050_d2-d1","V050_d8-d1","V050_d22-d1","V051_d8-d1","V051_d22-d1","V051_d2-d1","V054_d22-d1","V054_d2-d1","V054_d8-d1","V056_d2-d1","V056_d22-d1","V056_d8-d1","V057_d2-d1","V057_d22-d1","V057_d8-d1","V058_d2-d1","V058_d8-d1","V058_d22-d1","V061_d22-d1","V061_d2-d1","V061_d8-d1","V063_d22-d1","V063_d8-d1","V063_d2-d1","V064_d22-d1","V064_d8-d1","V064_d2-d1","V066_d8-d1","V066_d2-d1","V066_d22-d1","V067_d8-d1","V067_d22-d1","V067_d2-d1","V068_d22-d1","V068_d8-d1","V068_d2-d1","V070_d22-d1","V070_d8-d1","V070_d2-d1","V071_d8-d1","V071_d22-d1","V071_d2-d1","V073_d22-d1","V073_d2-d1","V073_d8-d1","V077_d22-d1","V077_d2-d1","V077_d8-d1","V140_d22-d1","V140_d8-d1","V140_d2-d1","V141_d8-d1","V141_d22-d1","V141_d2-d1","V142_d22-d1","V142_d2-d1","V142_d8-d1","V143_d2-d1","V143_d8-d1","V143_d22-d1","V145_d8-d1","V145_d22-d1","V145_d2-d1","V147_d8-d1","V147_d22-d1","V147_d2-d1","V148_d8-d1","V148_d22-d1","V148_d2-d1")


df.orig <- df %>% dplyr::select(one_of(c("Gene",samples_to_include)))
df.mat = df.orig[ , (colnames(df.orig) != "Gene" )] %>% as.data.frame
row.names(df.mat) <- df.orig$Gene
df.mat <- as.data.frame(df.mat)


annot <- read.csv("../results/naive_plot_metadata.csv", header = TRUE)
annot %>% dplyr::filter(sample_id %in% samples_to_include) -> annot
annot$contrast <- factor(annot$contrast)
groups = c("contrast")
relevel_factors <- TRUE
factor_relevel <- c("contrast:d2_d1,d8_d1,d22_d1,d23_d22")
if(relevel_factors){
  for(f in factor_relevel){
    variable <- unlist(str_split(f,":"))[1]
    if(!variable %in% groups){ next }
    levels <- unlist(str_split(unlist(str_split(f,":"))[2],","))
    annot[,variable] <- factor(as.character(annot[,variable]), levels = levels)
  }
}

if(TRUE){
  annot %>% dplyr::arrange_(.dots=groups) -> annot
  df.mat <- df.mat[,match(annot$sample_id,colnames(df.mat))] 
}

annot %>% dplyr::select(groups) -> annotation_col 
annotation_col = as.data.frame(unclass(annotation_col))
annotation_col[] <- lapply(annotation_col,factor)
rownames(annotation_col) <- annot$sample_id
annot_col = list()
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
qual_col_pals = qual_col_pals[c(7,6,2,1,8,3,4,5),]
colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

sample_color_palette <- FALSE
set_color_seed = FALSE
if(sample_color_palette || set_color_seed){
  seed <- if(set_color_seed){ 
    1234
  }else{
    sample(1:2^15,1)
  }
  
  cat(paste0("Using color seed: ",seed,"\n"))
  set.seed(seed)
  colors <- sample(colors)
}
#cat("Annotation color palette:\n")
#print(colors)
b=1
i=1
while (i <= length(groups)){
  nam <- groups[i]
  grp <- as.factor(annotation_col[,i])
  c <- b+length(levels(grp))-1
  col = colors[b:c]
  names(col) <- levels(grp)
  assign(nam,col)
  annot_col = append(annot_col,mget(nam))
  b = b+c
  i=i+1
}
cat("Annotation color palette:\n")
print(annot_col)

print(paste0("The total number of genes in heatmap: ", nrow(df.mat)))

labels_col <- colnames(df.mat)

imageWidth = 3000
imageHeight = 1500
dpi = 300

png(
  filename="../plots/naive_heatmap_vaccination_I.png",
  width=imageWidth,
  height=imageHeight,
  units="px",
  pointsize=4,
  bg="white",
  res=dpi,
  type="cairo")

p = doheatmap(dat=df.mat, clus=FALSE, clus2=TRUE, ht=50, rn=TRUE, cn=FALSE, col="Blue-Red 3",scale_range=1)

dev.off()
