# Functions defined here will be available to call in
# the code for any table.

install_bioconductor_package <- function(pkg) {
    install.packages(paste0("https://gypsum.palantircloud.com/assets/dyn/bioconductor-packages/", pkg, ".tar.gz"), repos=NULL)
}

@transform_pandas(
    Output(rid="ri.vector.main.execute.fe0728dc-02d9-4489-8712-3ff7e745ecb9"),
    combined_complex_metadata=Input(rid="ri.foundry.main.dataset.8b5c5b80-01ad-4311-9156-3aad6c944a48"),
    merged_diff=Input(rid="ri.foundry.main.dataset.7ea602f6-f92c-4696-8f77-b5d39152978e")
)
Vaccination_I_filt_heatmap <- function(merged_diff,combined_complex_metadata) {
    #This function uses pheatmap to draw a heatmap, scaling first by rows
    #(with samples in columns and genes in rows)
    # image: png
    suppressMessages(library(colorspace))
    suppressMessages(library(dendsort))
    suppressMessages(library(pheatmap))
    suppressMessages(library(tidyverse))
    suppressMessages(library(RColorBrewer))
    suppressMessages(library(dplyr))

    doheatmap <- function(dat, clus, clus2, ht, rn, cn, col) {
        require(pheatmap)
        require(dendsort)

        if (F) {
            tmean.scale = t(scale(t(dat)))
            tmean.scale = tmean.scale[!is.infinite(rowSums(tmean.scale)),]
            tmean.scale = na.omit(tmean.scale)
        } else {
            tmean.scale=dat
            #tmean.scale = na.omit(tmean.scale)
        }

        col_filter <- apply(as.data.frame(tmean.scale), 2, function(z) all(is.na(z)))
        row_filter <- apply(as.data.frame(tmean.scale), 1, function(z) length(z[is.na(z)]) > 10)
        tmean.scale <- tmean.scale[!row_filter, !col_filter]

        col.pal <- diverging_hcl(n=100, palette=col)
        if (FALSE) {
            col.pal = rev(col.pal)
        }
        # define metrics for clustering
        drows1 <- "correlation"
        dcols1 <- "correlation"

        minx = floor(min(tmean.scale))
        maxx = ceiling(max(tmean.scale))

        if (FALSE) {
            if(maxx > abs(minx)){
                minx = -1 * maxx
            }else{
                maxx = -1 * minx
            }
            breaks = seq(minx, maxx, length=100)
            legbreaks = seq(minx, maxx, length=5)
        } else {
            breaks = seq(-1, 1, length=100)
            legbreaks = seq(-1, 1, length=5)
        }
        breaks = sapply(breaks, signif, 4)
        legbreaks = sapply(legbreaks, signif, 4)

        #Run cluster method using 
        # hc = hclust(dist(t(tmean.scale)), method="average")
        # hcrow = hclust(dist(tmean.scale), method="average")

        # if (FALSE) {
        #     sort_hclust <- function(...) as.hclust(rev(dendsort(as.dendrogram(...))))
        # } else {
        #     sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
        # }
        # if (clus) {
        #     colclus <- sort_hclust(hc)
        # } else {
        #     colclus = FALSE
        # }
        # if (clus2) {
        #     rowclus <- sort_hclust(hcrow)
        # } else {
        #     rowclus = FALSE
        # }
        #print('sorted the clusters')
        if (TRUE) {
            treeheight <- 25
        } else {
            treeheight <- 0
        }

        gaps_col <- ""
        gaps_row <- ""

        gaps_col <- as.numeric(unlist(str_split(gaps_col,",")))
        gaps_row <- as.numeric(unlist(str_split(gaps_row,",")))

        gaps_col[is.na(gaps_col)] <- 0
        gaps_row[is.na(gaps_row)] <- 0

        hm.parameters <- list(
            tmean.scale, 
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
            gaps_col = gaps_col,
            gaps_row = gaps_row,
            annotation_names_col = FALSE,
            na_col = "#000000"
        )

        mat = t(tmean.scale)
        # print('calculated mat')
        
        do.call("pheatmap", c(hm.parameters))
    }

    df <- merged_diff

    samples_to_include = c("V002_d2-d1","V002_d22-d1","V002_d8-d1","V004_d22-d1","V004_d8-d1","V004_d2-d1","V005_d22-d1","V005_d8-d1","V005_d2-d1","V006_d22-d1","V006_d8-d1","V006_d2-d1","V007_d8-d1","V007_d2-d1","V007_d22-d1","V008_d8-d1","V008_d2-d1","V008_d22-d1","V009_d22-d1","V009_d8-d1","V009_d2-d1","V010_d22-d1","V010_d8-d1","V010_d2-d1","V011_d8-d1","V011_d22-d1","V011_d2-d1","V013_d22-d1","V013_d2-d1","V013_d8-d1","V015_d22-d1","V015_d2-d1","V015_d8-d1","V016_d22-d1","V016_d8-d1","V016_d2-d1","V017_d8-d1","V017_d22-d1","V017_d2-d1","V019_d22-d1","V019_d2-d1","V019_d8-d1","V020_d2-d1","V020_d8-d1","V020_d22-d1","V021_d8-d1","V021_d22-d1","V021_d2-d1","V022_d8-d1","V022_d22-d1","V022_d2-d1","V024_d8-d1","V024_d22-d1","V024_d2-d1","V025_d8-d1","V025_d22-d1","V025_d2-d1","V026_d8-d1","V026_d22-d1","V026_d2-d1","V027_d8-d1","V027_d22-d1","V027_d2-d1","V028_d8-d1","V028_d22-d1","V028_d2-d1","V029_d2-d1","V029_d8-d1","V029_d22-d1","V030_d8-d1","V030_d22-d1","V030_d2-d1","V031_d22-d1","V031_d8-d1","V031_d2-d1","V033_d22-d1","V033_d8-d1","V033_d2-d1","V034_d8-d1","V034_d22-d1","V034_d2-d1","V036_d2-d1","V036_d8-d1","V036_d22-d1","V037_d22-d1","V037_d8-d1","V037_d2-d1","V039_d22-d1","V039_d8-d1","V039_d2-d1","V041_d22-d1","V041_d2-d1","V041_d8-d1","V043_d2-d1","V043_d8-d1","V043_d22-d1","V045_d22-d1","V045_d2-d1","V045_d8-d1","V048_d22-d1","V048_d2-d1","V048_d8-d1","V049_d22-d1","V049_d2-d1","V049_d8-d1","V050_d2-d1","V050_d8-d1","V050_d22-d1","V051_d8-d1","V051_d22-d1","V051_d2-d1","V054_d22-d1","V054_d2-d1","V054_d8-d1","V056_d2-d1","V056_d22-d1","V056_d8-d1","V057_d2-d1","V057_d22-d1","V057_d8-d1","V058_d2-d1","V058_d8-d1","V058_d22-d1","V061_d22-d1","V061_d2-d1","V061_d8-d1","V063_d22-d1","V063_d8-d1","V063_d2-d1","V064_d22-d1","V064_d8-d1","V064_d2-d1","V066_d8-d1","V066_d2-d1","V066_d22-d1","V067_d8-d1","V067_d22-d1","V067_d2-d1","V068_d22-d1","V068_d8-d1","V068_d2-d1","V070_d22-d1","V070_d8-d1","V070_d2-d1","V071_d8-d1","V071_d22-d1","V071_d2-d1","V073_d22-d1","V073_d2-d1","V073_d8-d1","V077_d22-d1","V077_d2-d1","V077_d8-d1","V140_d22-d1","V140_d8-d1","V140_d2-d1","V141_d8-d1","V141_d22-d1","V141_d2-d1","V142_d22-d1","V142_d2-d1","V142_d8-d1","V143_d2-d1","V143_d8-d1","V143_d22-d1","V145_d8-d1","V145_d22-d1","V145_d2-d1","V147_d8-d1","V147_d22-d1","V147_d2-d1","V148_d8-d1","V148_d22-d1","V148_d2-d1")

    
    df.orig <- df %>% dplyr::select(one_of(c("Gene",samples_to_include)))
    df.mat = df.orig[ , (colnames(df.orig) != "Gene" )] %>% as.data.frame
    row.names(df.mat) <- df.orig$Gene
    df.mat <- as.data.frame(df.mat)
    
    
    annot <- combined_complex_metadata
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
      filename=graphicsFile,
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")

    p = doheatmap(dat=df.mat, clus=FALSE, clus2=TRUE, ht=50, rn=TRUE, cn=FALSE, col="Blue-Red 3")
    #Order output    
    row.order <- if("order" %in% names(p$tree_row)){
        p$tree_row[["order"]]
    }else{
        row.names(df.mat)
    }

    col.order <- if("order" %in% names(p$tree_col)){
        p$tree_col[["order"]]
    }else{
        colnames(df.mat)
    }

    df.mat <- df.mat[row.order,col.order]
    df.mat <- df.mat %>% rownames_to_column("Gene")

    return(df.mat)
}

@transform_pandas(
    Output(rid="ri.vector.main.execute.497dca58-690d-4b21-824e-63ea4ad797e9"),
    combined_complex_metadata_immune=Input(rid="ri.foundry.main.dataset.7318826e-2f14-4395-8720-f2af0556bef7"),
    merged_diff_subset_immune=Input(rid="ri.foundry.main.dataset.a702eeb4-44cf-4dc4-85a5-462e19c6563f")
)
Vaccination_I_filt_heatmap_copied <- function(merged_diff_subset_immune,combined_complex_metadata_immune) {
    #This function uses pheatmap to draw a heatmap, scaling first by rows
    #(with samples in columns and genes in rows)
    # image: png
    suppressMessages(library(colorspace))
    suppressMessages(library(dendsort))
    suppressMessages(library(pheatmap))
    suppressMessages(library(tidyverse))
    suppressMessages(library(RColorBrewer))
    suppressMessages(library(dplyr))

    doheatmap <- function(dat, clus, clus2, ht, rn, cn, col) {
        require(pheatmap)
        require(dendsort)

        if (F) {
            tmean.scale = t(scale(t(dat)))
            tmean.scale = tmean.scale[!is.infinite(rowSums(tmean.scale)),]
            tmean.scale = na.omit(tmean.scale)
        } else {
            tmean.scale=dat
            #tmean.scale = na.omit(tmean.scale)
        }
        col_filter <- apply(as.data.frame(tmean.scale), 2, function(z) all(is.na(z)))
        row_filter <- apply(as.data.frame(tmean.scale), 1, function(z) length(z[is.na(z)]) > 5)
        tmean.scale <- tmean.scale[!row_filter, !col_filter]

        col.pal <- diverging_hcl(n=100, palette=col)
        if (FALSE) {
            col.pal = rev(col.pal)
        }
        # define metrics for clustering
        drows1 <- "correlation"
        dcols1 <- "correlation"

        minx = floor(min(tmean.scale))
        maxx = ceiling(max(tmean.scale))

        if (FALSE) {
            if(maxx > abs(minx)){
                minx = -1 * maxx
            }else{
                maxx = -1 * minx
            }
            breaks = seq(minx, maxx, length=100)
            legbreaks = seq(minx, maxx, length=5)
        } else {
            breaks = seq(-2, 2, length=100)
            legbreaks = seq(-2, 2, length=5)
        }
        breaks = sapply(breaks, signif, 4)
        legbreaks = sapply(legbreaks, signif, 4)

        if (TRUE) {
            treeheight <- 25
        } else {
            treeheight <- 0
        }

        gaps_col <- ""
        gaps_row <- ""

        gaps_col <- as.numeric(unlist(str_split(gaps_col,",")))
        gaps_row <- as.numeric(unlist(str_split(gaps_row,",")))

        gaps_col[is.na(gaps_col)] <- 0
        gaps_row[is.na(gaps_row)] <- 0

        hm.parameters <- list(
            tmean.scale, 
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
            gaps_col = gaps_col,
            gaps_row = gaps_row,
            annotation_names_col = FALSE,
            na_col = "#000000"
        )

        mat = t(tmean.scale)
        # print('calculated mat')
        
        do.call("pheatmap", c(hm.parameters))
    }

    df <- merged_diff_subset_immune

    samples_to_include = c("V001_d2-d1","V003_d2-d1","V014_d2-d1","V032_d2-d1","V062_d2-d1","V001_d8-d1","V003_d8-d1","V014_d8-d1","V032_d8-d1","V062_d8-d1","V001_d22-d1","V003_d22-d1","V014_d22-d1","V032_d22-d1","V062_d22-d1")

    df.orig <- df %>% dplyr::select(one_of(c("Gene",samples_to_include)))
    df.mat = df.orig[ , (colnames(df.orig) != "Gene" )] %>% as.data.frame
    row.names(df.mat) <- df.orig$Gene
    df.mat <- as.data.frame(df.mat)
    
    annot <- combined_complex_metadata_immune
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
      filename=graphicsFile,
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")

    p = doheatmap(dat=df.mat, clus=FALSE, clus2=TRUE, ht=50, rn=TRUE, cn=FALSE, col="Blue-Red 3")
    #Order output    
    row.order <- if("order" %in% names(p$tree_row)){
        p$tree_row[["order"]]
    }else{
        row.names(df.mat)
    }

    col.order <- if("order" %in% names(p$tree_col)){
        p$tree_col[["order"]]
    }else{
        colnames(df.mat)
    }

    df.mat <- df.mat[row.order,col.order]
    df.mat <- df.mat %>% rownames_to_column("Gene")

    return(df.mat)
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.8b5c5b80-01ad-4311-9156-3aad6c944a48"),
    diff_diff_metadata=Input(rid="ri.foundry.main.dataset.30626dd8-9bd4-444c-b1c9-48c4de1e9e83"),
    metadata=Input(rid="ri.foundry.main.dataset.76aa3c8e-4e1d-4978-966b-41633299a883")
)
combined_complex_metadata <- function(diff_diff_metadata, metadata) {

    complex.df <- diff_diff_metadata
    merged.df <- metadata %>% select(-sample_id,-timepoint) %>% distinct()
    print(head(merged.df))
    df <- merge(complex.df,merged.df, by="patient_id")
    df <- df %>% select(sample_id, everything()) %>% distinct() 

    return(df)   
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.7318826e-2f14-4395-8720-f2af0556bef7"),
    complex_metadata_copied=Input(rid="ri.foundry.main.dataset.c9121cc9-ae63-496c-869f-adc0039d507e"),
    metadata=Input(rid="ri.foundry.main.dataset.76aa3c8e-4e1d-4978-966b-41633299a883")
)
combined_complex_metadata_immune <- function(complex_metadata_copied, metadata) {

    complex.df <- complex_metadata_copied
    merged.df <- metadata %>% select(-sample_id,-timepoint) %>% distinct()
    print(head(merged.df))
    df <- merge(complex.df,merged.df, by="patient_id")
    df <- df %>% select(sample_id, everything()) %>% distinct() 

    return(df)   
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.2bdc84b9-8014-4e76-bf11-19e4dcef5968"),
    naive_diff_counts=Input(rid="ri.foundry.main.dataset.a1fcac5c-878d-4d65-a76a-bd9c9794fb10"),
    naive_diff_metadata=Input(rid="ri.foundry.main.dataset.2e114f09-86e4-4d4f-9684-44ac9e05d4e1")
)
complex_diff <- function(naive_diff_counts, naive_diff_metadata) {

    library(tidyverse)

    df.orig <- naive_diff_counts
    annot <- naive_diff_metadata
    
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
    
    df.ret <- df.ret %>% rownames_to_column("Gene")
    return(df.ret)
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.c01d46dc-9e19-4903-8861-8aa3d5b7f8a2"),
    fc_metadata_immune=Input(rid="ri.foundry.main.dataset.7171519d-c894-4431-8046-83701c07fa5b"),
    grouped_counts_immune=Input(rid="ri.foundry.main.dataset.f81a9315-a47d-4f42-8d1f-cce65b2e920f")
)
complex_diff_immune <- function(grouped_counts_immune, fc_metadata_immune) {

    library(tidyverse)

    df.orig <- grouped_counts_immune
    annot <- fc_metadata_immune
    
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
    
    df.ret <- df.ret %>% rownames_to_column("Gene")
    return(df.ret)
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.c9121cc9-ae63-496c-869f-adc0039d507e"),
    complex_diff_immune=Input(rid="ri.foundry.main.dataset.c01d46dc-9e19-4903-8861-8aa3d5b7f8a2"),
    fc_metadata_immune=Input(rid="ri.foundry.main.dataset.7171519d-c894-4431-8046-83701c07fa5b")
)
complex_metadata_copied <- function(complex_diff_immune, fc_metadata_immune) {
    
    library(stringr)

    df <- complex_diff_immune
    annot <- fc_metadata_immune

    col.nam <- colnames(df)[colnames(df) != "Gene"]
    
    patient_id <- apply(array(col.nam),1,function(z) unlist(str_split(z,"_"))[1]) #  unlist(str_split(col.nam,"_"))[1]
    contrast <- rep("vac2_vac1",length(patient_id)) #apply(array(col.nam),1,function(z) unlist(str_split(z,"_"))[2])
    #t.level <- apply(array(contrast),1,function(z) unlist(str_split(z,"_"))[1]) #unlist(str_split(contrast,"-"))[1]
    #r.level <- apply(array(contrast),1,function(z) unlist(str_split(z,"_"))[2]) #unlist(str_split(contrast,"-"))[2]

    df.ret <- data.frame(sample_id=col.nam, patient_id = patient_id, contrast = gsub("-","_",contrast))
    df <- rbind(annot,df.ret)
    return(df)
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.107ded67-aa74-4f7b-b01f-a950bc011336"),
    merged_diff=Input(rid="ri.foundry.main.dataset.7ea602f6-f92c-4696-8f77-b5d39152978e")
)
corr_plot_d23vd22 <- function(merged_diff) {
    #image: png
    imageWidth = 1300
    imageHeight = 1300
    dpi = 300

    png(
      filename=graphicsFile,
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")
   
    library(corrplot)
    library(tidyverse)
    library(colorspace)

    df <- merged_diff

    contrast <- "d23-d22$"

    sample_of_interest <- grep(contrast, colnames(df), value=TRUE)

    df <- df %>% select(Gene,one_of(sample_of_interest))
    row.names(df) <- df$Gene
    df$Gene <- NULL

    converter_in  <- c("alpha","beta","gamma","delta","IL-12IL-23p40")
    converter_out <- c("a","b","g","d","p40")
    
    rownames(df) <- .multigsub(converter_in,converter_out,rownames(df))

    df <- as.data.frame(t(as.matrix(df)))
    #remove bad columns
    filter <- apply(df, 2, function(z) all(is.na(z)))
    df <- df[,!filter]
    # Remove
    filter <- apply(df,1,function(z) any(is.na(z)))
    df <- df[!filter,]

    M <- cor(df, method = "spearman")
    p.mat <- cor.mtest(df)$p
    p.mat <- matrix(p.adjust(c(p.mat)), nrow=nrow(p.mat), byrow=FALSE)
    head(p.mat[, 1:5])
    col <- diverging_hcl(n=200, palette="Blue-Red")
    
    p <- corrplot(M, method = "ellipse", col = col,
         type = "upper", order = "hclust", number.cex = 1, number.font=2,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, tl.cex=1.5, # Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag = FALSE,cl.cex = 2)
    # p <- .corrplot.mixed(M, #col = col(200),
    #      upper="color", lower="ellipse", order = "hclust", number.cex = 2, number.font=2, diag = "n",
    #      addCoef.col = "black", # Add coefficient of correlation
    #      tl.col = "black", tl.srt = 90, tl.cex=2, # Text label color and rotation
    #      # Combine with significance
    #      p.mat = p.mat, sig.level = 0.05, insig = "blank")

    #corrplot.mixed(M,upper="color",lower="ellipse")
    # p <- corrplot.mixed(M, #col = col(200),
    #      upper="color", lower="ellipse", order = "hclust", number.cex = 2, number.font=2, diag = "n",
    #      addCoef.col = "black", # Add coefficient of correlation
    #      tl.col = "black", tl.srt = 90, tl.cex=2, # Text label color and rotation
    #      # Combine with significance
    #      p.mat = p.mat, sig.level = 0.05, insig = "blank", upper.col=col,lower.col=col, cl.cex = 3)
    df <- rownames_to_column(df, "patient_id")
    return(df)
}

.multigsub <-
function (pattern, replacement, text.var, leadspace = FALSE, 
    trailspace = FALSE, fixed = TRUE, trim = TRUE, order.pattern = fixed, 
    ...) {

    if (leadspace | trailspace) replacement <- spaste(replacement, trailing = trailspace, leading = leadspace)

    if (fixed && order.pattern) {
        ord <- rev(order(nchar(pattern)))
        pattern <- pattern[ord]
        if (length(replacement) != 1) replacement <- replacement[ord]
    }
    if (length(replacement) == 1) replacement <- rep(replacement, length(pattern))
   
    for (i in seq_along(pattern)){
        text.var <- gsub(pattern[i], replacement[i], text.var, fixed = fixed, ...)
    }

    if (trim) text.var <- gsub("\\s+", " ", gsub("^\\s+|\\s+$", "", text.var, perl=TRUE), perl=TRUE)
    text.var
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.196fe3c5-0f03-4e47-a161-52078284c41b"),
    merged_diff=Input(rid="ri.foundry.main.dataset.7ea602f6-f92c-4696-8f77-b5d39152978e")
)
corr_plot_d2vd1 <- function(merged_diff) {
    #image: png
    imageWidth = 1300
    imageHeight = 1300
    dpi = 300

    png(
      filename=graphicsFile,
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")
              
    library(corrplot)
    library(tidyverse)
    library(colorspace)

    df <- merged_diff

    contrast <- "_d2-d1$"

    sample_of_interest <- grep(contrast, colnames(df), value=TRUE)

    df <- df %>% select(Gene,one_of(sample_of_interest))
    row.names(df) <- df$Gene
    df$Gene <- NULL

    converter_in  <- c("alpha","beta","gamma","delta","IL-12IL-23p40")
    converter_out <- c("a","b","g","d","p40")
    
    rownames(df) <- .multigsub(converter_in,converter_out,rownames(df))
    
    df <- as.data.frame(t(as.matrix(df)))

    # Remove
    #remove bad columns
    filter <- apply(df, 2, function(z) all(is.na(z)))
    df <- df[,!filter]
    # Remove
    filter <- apply(df,1,function(z) any(is.na(z)))
    df <- df[!filter,]
    
    M <- cor(df, method = "spearman")
    p.mat <- cor.mtest(df)$p

    p.mat <- matrix(p.adjust(c(p.mat)), nrow=nrow(p.mat), byrow=FALSE)
    col <- diverging_hcl(n=200, palette="Blue-Red")

    p <- corrplot(M, method = "ellipse", col = col,
         type = "upper", order = "hclust", number.cex = 1, number.font=2,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, tl.cex=1.5, # Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag = FALSE,cl.cex = 2)
    # p <- .corrplot.mixed(M, #col = col(200),
    #      upper="color", lower="ellipse", order = "hclust", number.cex = 2, number.font=2, diag = "n",
    #      addCoef.col = "black", # Add coefficient of correlation
    #      tl.col = "black", tl.srt = 90, tl.cex=2, # Text label color and rotation
    #      # Combine with significance
    #      p.mat = p.mat, sig.level = 0.05, insig = "blank")

    #corrplot.mixed(M,upper="color",lower="ellipse")
    # p <- corrplot.mixed(M, #col = col(200),
    #      upper="color", lower="ellipse", order = "hclust", number.cex = 2, number.font=2, diag = "n",
    #      addCoef.col = "black", # Add coefficient of correlation
    #      tl.col = "black", tl.srt = 90, tl.cex=2, # Text label color and rotation
    #      # Combine with significance
    #      p.mat = p.mat, sig.level = 0.05, insig = "blank", upper.col=col,lower.col=col, cl.cex = 3)
    df <- rownames_to_column(df, "patient_id")
    return(df)
}

.multigsub <-
function (pattern, replacement, text.var, leadspace = FALSE, 
    trailspace = FALSE, fixed = TRUE, trim = TRUE, order.pattern = fixed, 
    ...) {

    if (leadspace | trailspace) replacement <- spaste(replacement, trailing = trailspace, leading = leadspace)

    if (fixed && order.pattern) {
        ord <- rev(order(nchar(pattern)))
        pattern <- pattern[ord]
        if (length(replacement) != 1) replacement <- replacement[ord]
    }
    if (length(replacement) == 1) replacement <- rep(replacement, length(pattern))
   
    for (i in seq_along(pattern)){
        text.var <- gsub(pattern[i], replacement[i], text.var, fixed = fixed, ...)
    }

    if (trim) text.var <- gsub("\\s+", " ", gsub("^\\s+|\\s+$", "", text.var, perl=TRUE), perl=TRUE)
    text.var
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.30626dd8-9bd4-444c-b1c9-48c4de1e9e83"),
    complex_diff=Input(rid="ri.foundry.main.dataset.2bdc84b9-8014-4e76-bf11-19e4dcef5968"),
    naive_diff_metadata=Input(rid="ri.foundry.main.dataset.2e114f09-86e4-4d4f-9684-44ac9e05d4e1")
)
diff_diff_metadata <- function(complex_diff, naive_diff_metadata) {
    
    library(stringr)

    df <- complex_diff
    annot <- naive_diff_metadata

    col.nam <- colnames(df)[colnames(df) != "Gene"]
    patient_id <- apply(array(col.nam),1,function(z) unlist(str_split(z,"_"))[1]) #  unlist(str_split(col.nam,"_"))[1]
    contrast <- rep("vac2_vac1",length(patient_id)) #apply(array(col.nam),1,function(z) unlist(str_split(z,"_"))[2])
    #t.level <- apply(array(contrast),1,function(z) unlist(str_split(z,"_"))[1]) #unlist(str_split(contrast,"-"))[1]
    #r.level <- apply(array(contrast),1,function(z) unlist(str_split(z,"_"))[2]) #unlist(str_split(contrast,"-"))[2]

    df.ret <- data.frame(sample_id=col.nam, patient_id = patient_id, contrast = gsub("-","_",contrast))
    df <- rbind(annot,df.ret)
    return(df)
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.7171519d-c894-4431-8046-83701c07fa5b"),
    grouped_counts_immune=Input(rid="ri.foundry.main.dataset.f81a9315-a47d-4f42-8d1f-cce65b2e920f")
)
fc_metadata_immune <- function(grouped_counts_immune) {
    library(stringr)
    df <- grouped_counts_immune

    col.nam <- colnames(df)[colnames(df) != "Gene"]
    
    patient_id <- apply(array(col.nam),1,function(z) unlist(str_split(z,"_"))[1]) #  unlist(str_split(col.nam,"_"))[1]
    contrast <- apply(array(col.nam),1,function(z) unlist(str_split(z,"_"))[2])
    #t.level <- apply(array(contrast),1,function(z) unlist(str_split(z,"_"))[1]) #unlist(str_split(contrast,"-"))[1]
    #r.level <- apply(array(contrast),1,function(z) unlist(str_split(z,"_"))[2]) #unlist(str_split(contrast,"-"))[2]

    df.ret <- data.frame(sample_id=col.nam, patient_id = patient_id, contrast = gsub("-","_",contrast))

    return(df.ret)

}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.f81a9315-a47d-4f42-8d1f-cce65b2e920f"),
    immune_dataset=Input(rid="ri.foundry.main.dataset.f932b974-708c-42f9-92da-bf1534039543"),
    metadata=Input(rid="ri.foundry.main.dataset.76aa3c8e-4e1d-4978-966b-41633299a883")
)
grouped_counts_immune <- function(immune_dataset, metadata) {
    
    library(tidyverse)
    df <- immune_dataset
    annot <- metadata

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
            df.all <- df.ret
        }else{
            df.all <- merge(df.all,df.ret,by="Gene")
        }
    }

    return(df.all)
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.f932b974-708c-42f9-92da-bf1534039543"),
    preimmune_patients=Input(rid="ri.foundry.main.dataset.3874744d-f071-47bb-9dcd-7f5dc80a8efa")
)
immune_dataset <- function(msd_dataset, preimmune_patients) {

    df <- msd_dataset
    preimmune <- preimmune_patients

    filter <- grep(paste(preimmune$patient_id,collapse="|"),colnames(df),value=TRUE)
    df <- df%>% select(one_of(c("Gene",filter)))
    return(df)
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.7ea602f6-f92c-4696-8f77-b5d39152978e"),
    complex_diff=Input(rid="ri.foundry.main.dataset.2bdc84b9-8014-4e76-bf11-19e4dcef5968"),
    naive_diff_counts=Input(rid="ri.foundry.main.dataset.a1fcac5c-878d-4d65-a76a-bd9c9794fb10")
)
merged_diff <- function(naive_diff_counts,complex_diff) {
    
    df1 <- naive_diff_counts
    df2 <- complex_diff

    df <- merge(df1,df2,by="Gene")
    
    return(df)
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.9123c7bd-0782-4c95-a454-4c297fbac4bb"),
    complex_diff_immune=Input(rid="ri.foundry.main.dataset.c01d46dc-9e19-4903-8861-8aa3d5b7f8a2"),
    grouped_counts_immune=Input(rid="ri.foundry.main.dataset.f81a9315-a47d-4f42-8d1f-cce65b2e920f")
)
merged_diff_immune <- function(grouped_counts_immune,complex_diff_immune) {
    
    df1 <- grouped_counts_immune
    df2 <- complex_diff_immune

    df <- merge(df1,df2,by="Gene")
    
    return(df)
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.a702eeb4-44cf-4dc4-85a5-462e19c6563f"),
    merged_diff_immune=Input(rid="ri.foundry.main.dataset.9123c7bd-0782-4c95-a454-4c297fbac4bb"),
    selected_analytes=Input(rid="ri.foundry.main.dataset.4ba63091-d7c1-4e82-b2b9-6fc18808feab")
)
merged_diff_subset_immune <- function(merged_diff_immune, selected_analytes) {
    
    analytes <- selected_analytes$Gene

    df <- merged_diff_immune

    print(analytes[!analytes %in% df$Gene])

    df <- df %>% filter(Gene %in% analytes)

    return(df)
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.78426f81-b289-4600-b6ea-d97ccf829f1d"),
    msd_dataset=Input(rid="ri.foundry.main.dataset.6004d02f-02c1-419f-ab0c-b19ebcd8c9aa"),
    preimmune_patients=Input(rid="ri.foundry.main.dataset.3874744d-f071-47bb-9dcd-7f5dc80a8efa")
)
naive_dataset <- function(msd_dataset, preimmune_patients) {

    df <- msd_dataset
    preimmune <- preimmune_patients

    filter <- grep(paste(preimmune$patient_id,collapse="|"),colnames(df),value=TRUE)
    df <- df%>% select(-one_of(filter))
    return(df)
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.862f8327-dcb8-414e-b8ae-fb74aaa3a15b"),
    metadata=Input(rid="ri.foundry.main.dataset.76aa3c8e-4e1d-4978-966b-41633299a883"),
    naive_dataset=Input(rid="ri.foundry.main.dataset.78426f81-b289-4600-b6ea-d97ccf829f1d")
)
naive_deg_analysis <- function(naive_dataset, metadata) {
    suppressMessages(library(limma))
    suppressMessages(library(tidyverse))
    suppressMessages(library(edgeR))
    suppressMessages(library(stringr))
    
    df <- naive_dataset
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
    
    targetfile <- metadata
    targetfile <- targetfile[match(colnames(df.m),targetfile$sample_id),]
    targetfile <- targetfile[rowSums(is.na(targetfile)) != ncol(targetfile), ]
    df.m <- df.m[,match(targetfile$sample_id,colnames(df.m))]

    x <- 2^df.m

    ordered_covariates=c("timepoint","patient_id")
    
    ordered_covariates=ordered_covariates[order(ordered_covariates!="timepoint")]

    print(tail(targetfile))
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
                                "d22-d1",
                                "d8-d1",
                                "(d23-d22)-(d2-d1)")

    cm <- makeContrasts(contrasts = contrasts_of_interest, levels=design)
    
    fit2 <- contrasts.fit(fit, cm)

    fit2 <- eBayes(fit2)

    logFC = fit2$coefficients
    colnames(logFC)=paste(colnames(logFC),"logFC",sep="_")
    tstat = fit2$t
    colnames(tstat)=paste(colnames(tstat),"tstat",sep="_")
    FC = 2^fit2$coefficients
    FC = ifelse(FC<1,-1/FC,FC)
    colnames(FC)=paste(colnames(FC),"FC",sep="_")
    pvalall=fit2$p.value
    colnames(pvalall)=paste(colnames(pvalall),"pval",sep="_")
    pvaladjall=apply(pvalall,2,function(x) p.adjust(x,"BH"))
    colnames(pvaladjall)=paste(colnames(fit2$coefficients),"adjpval",sep="_")
    
    finalres=as.data.frame(cbind(v$E,FC, logFC, tstat, pvalall, pvaladjall))

    finalres %>% rownames_to_column("Gene") -> finalres
    print(paste0("Total number of genes included: ", nrow(finalres)))
    
    call_me_alias<-colnames(finalres)
    colnames(finalres)<-gsub("\\(|\\)","",call_me_alias)
    sparkdffr<-createDataFrame(finalres)
    return(sparkdffr) 
}

#################################################
## Global imports and functions included below ##
#################################################

# voom function modified to remove lib.size normalization and normalizeBetweenArrays steps
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

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.a1fcac5c-878d-4d65-a76a-bd9c9794fb10"),
    genes_of_interest=Input(rid="ri.foundry.main.dataset.a281bca8-cda6-4d1f-b03a-b98199d9b7b1"),
    metadata=Input(rid="ri.foundry.main.dataset.76aa3c8e-4e1d-4978-966b-41633299a883"),
    naive_dataset=Input(rid="ri.foundry.main.dataset.78426f81-b289-4600-b6ea-d97ccf829f1d")
)
naive_diff_counts <- function(naive_dataset, metadata, genes_of_interest) {
    
    library(tidyverse)
    df <- naive_dataset
    annot <- metadata
    goi <- genes_of_interest

    df <- df %>% filter(Gene %in% goi$Gene)
    
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
            df.all <- df.ret
        }else{
            df.all <- merge(df.all,df.ret,by="Gene")
        }
    }

    return(df.all)
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.2e114f09-86e4-4d4f-9684-44ac9e05d4e1"),
    naive_diff_counts=Input(rid="ri.foundry.main.dataset.a1fcac5c-878d-4d65-a76a-bd9c9794fb10")
)
naive_diff_metadata <- function(naive_diff_counts) {
    library(stringr)
    df <- naive_diff_counts

    col.nam <- colnames(df)[colnames(df) != "Gene"]
    
    patient_id <- apply(array(col.nam),1,function(z) unlist(str_split(z,"_"))[1]) #  unlist(str_split(col.nam,"_"))[1]
    contrast <- apply(array(col.nam),1,function(z) unlist(str_split(z,"_"))[2])
    #t.level <- apply(array(contrast),1,function(z) unlist(str_split(z,"_"))[1]) #unlist(str_split(contrast,"-"))[1]
    #r.level <- apply(array(contrast),1,function(z) unlist(str_split(z,"_"))[2]) #unlist(str_split(contrast,"-"))[2]

    df.ret <- data.frame(sample_id=col.nam, patient_id = patient_id, contrast = gsub("-","_",contrast))

    return(df.ret)

}

@transform_pandas(
    Output(rid="ri.vector.main.execute.ea7afa19-228e-4ec9-8776-1505197408f6"),
    corr_plot_d2vd1=Input(rid="ri.foundry.main.dataset.196fe3c5-0f03-4e47-a161-52078284c41b")
)
scatter_plot_v1_noline_skinny <- function(corr_plot_d2vd1) {
    #image: png
    imageWidth = (1300/4)*3
    imageHeight = (1300/4)*3
    dpi = 300

    png(
      filename=graphicsFile,
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")
          
    library(gridExtra)
    library(grid)
    library(ggplot2)
    library(lattice)
    library(stringr)
    library(ggpubr)

    df <- corr_plot_d2vd1 %>% select(-patient_id)

    colnames(df) <- gsub("-","", colnames(df))

    # IL6 (x) vs IP-10 (y)
    # TNFa (x) vs IP-10 (y)
    # TNFa (x) vs IL-6 (y)
    # IFNg (x) vs IP-10 (y)
    # IFNg (x) vs IL-6 (y)
    # IFNg (x) vs TNFa (y)
    # IL-15 (x) vs IP-10 (y)
    # IL-15 (x) vs IL-6 (y)
    # IL-15 (x) vs TNFa (y)
    # IL-15 (x) vs IFN (y)
    plots <- c("IL15-IFNg",
               "IL15-IL6", "IFNg-IL6",
               "IL15-IP10","IFNg-IP10", "IL6-IP10")
    lay <- rbind(c( 1, NA, NA),
                 c( 2, 3,  NA),
                 c( 4, 5,   6))

    plot_cor <- function(name){
        text.size = 6

        x <- unlist(str_split(name,"-"))[1]
        y <- unlist(str_split(name,"-"))[2]

        if(x == "IL15"){
            axis.title.y = element_text(size = text.size+2, colour="black")
            axis.text.y=element_text(colour = "black", size=text.size)
        }else{
            axis.title.y = element_text(size = 2, colour="white")
            axis.text.y=element_text(colour = "white", size=1)
        }

        if(y == "IP10"){
            axis.title.x = element_text(size = text.size+2, colour="black")
            axis.text.x=element_text(colour="black", size=text.size)
        }else{
            axis.title.x = element_text(size = text.size+2, colour="white")
            axis.text.x=element_text(colour="white", size=text.size)
        }

        plt <- ggscatter(df, x=x, y=y, cor.method="spearman", color = "black", shape = 19, size = 0.2, # Points color, shape and size

                    conf.int = FALSE, # Add confidence interval
                    cor.coef = FALSE)+
                    theme(  axis.title.x = axis.title.x,
                            axis.title.y = axis.title.y,
                            axis.text.x = axis.text.x,
                            axis.text.y = axis.text.y,
                            axis.line = element_line(size=0.2),
                            axis.ticks= element_line(size=0.2))+
                    scale_y_continuous(labels = scales::number_format(accuracy = 0.1, trim=TRUE))+
                    scale_x_continuous(labels = scales::number_format(accuracy = 0.1, trim=TRUE))

        return(plt)
    }

    grob.list <- lapply(plots,plot_cor)
    select_grobs <- function(hlay) {
        id <- unique(c(t(hlay))) 
        id[!is.na(id)]
    } 
    p <- grid.arrange(grobs = grob.list[select_grobs(lay)], layout_matrix = lay)

    print(p)
    return(NULL)
}

@transform_pandas(
    Output(rid="ri.vector.main.execute.bb3efd3a-57db-400e-92ac-85b438f3437b"),
    corr_plot_d23vd22=Input(rid="ri.foundry.main.dataset.107ded67-aa74-4f7b-b01f-a950bc011336")
)
scatter_plot_v2_alt_noline_skinny <- function(corr_plot_d23vd22) {
    #image: png
    imageWidth = (1300/4)*6
    imageHeight = (1300/4)*3
    dpi = 300

    png(
      filename=graphicsFile,
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")
          
    library(gridExtra)
    library(grid)
    library(ggplot2)
    library(lattice)
    library(stringr)
    library(ggpubr)

    df <- corr_plot_d23vd22 %>% select(-patient_id)

    colnames(df) <- gsub("-","", colnames(df))

    # IL6 (x) vs IP-10 (y)
    # TNFa (x) vs IP-10 (y)
    # TNFa (x) vs IL-6 (y)
    # IFNg (x) vs IP-10 (y)
    # IFNg (x) vs IL-6 (y)
    # IFNg (x) vs TNFa (y)
    # IL-15 (x) vs IP-10 (y)
    # IL-15 (x) vs IL-6 (y)
    # IL-15 (x) vs TNFa (y)
    # IL-15 (x) vs IFN (y)
    plots <- c("MCP1-IL1Ra", "MIP1b-IL1Ra", "IL6-IL1Ra", "TNFa-IL1Ra", "IFNg-IL1Ra", "IL15-IL1Ra",
               "MCP1-MCP1",  "MIP1b-MCP1",  "IL6-MCP1",  "TNFa-MCP1",  "IFNg-MCP1",  "IL15-MCP1",
               "MCP1-MIP1b", "MIP1b-MIP1b", "IL6-MIP1b", "TNFa-MIP1b", "IFNg-MIP1b", "IL15-MIP1b")
    lay <- rbind(c( 1,  2,  3,  4,  5,  6),
                 c( 7,  8,  9, 10, 11, 12),
                 c(13, 14, 15, 16, 17, 18))

    plot_cor <- function(i){
        text.size = 4.2
        name <- plots[i]
        x <- unlist(str_split(name,"-"))[1]
        y <- unlist(str_split(name,"-"))[2]

        axis.text.x=element_text(colour="black", size=text.size)
        axis.text.y=element_text(colour = "black", size=text.size)

        if(x == "MCP1"){
            axis.title.y = element_text(size = text.size+2, colour="black")
        }else{
            axis.title.y = element_text(size = text.size+2, colour="white")
        }

        if(y == "MIP1b"){
            axis.title.x = element_text(size = text.size+2, colour="black")
        }else{
            axis.title.x = element_text(size = text.size+2, colour="white")
        }

        if(i %in% c(1:6,8,16,18)){
            plt <- ggscatter(df, x=x, y=y, cor.method="spearman", color = "black", shape = 19, size = 0.2, # Points color, shape and size

                        conf.int = FALSE, # Add confidence interval
                        cor.coef = FALSE)+
                    theme(  axis.title.x = axis.title.x,
                            axis.title.y = axis.title.y,
                            axis.text.x = axis.text.x,
                            axis.text.y = axis.text.y,
                            axis.line = element_line(size=0.2),
                            axis.ticks= element_line(size=0.2))
        }else{
            #Print blanks with labels
            axis.text = element_text(colour = "white", size=1)
            axis.line = element_line(colour = "white", size=0.2)
            axis.ticks = element_line(colour= "white", size=0.2)
            plt <- ggscatter(df, x=x, y=y, cor.method="spearman", color = "white", shape = 19, size = 0.2, # Points color, shape and size

                        conf.int = FALSE, # Add confidence interval
                        cor.coef = FALSE)+
                    theme(  axis.title.x = axis.title.x,
                            axis.title.y = axis.title.y,
                            axis.text = axis.text,
                            axis.line = axis.line,
                            axis.ticks= axis.ticks)           
        }

        plt <- plt + scale_y_continuous(labels = scales::number_format(accuracy = 0.1, trim=TRUE))+
                     scale_x_continuous(labels = scales::number_format(accuracy = 0.1, trim=TRUE))

        return(plt)
    }

    grob.list <- lapply(seq_along(plots),plot_cor)
    select_grobs <- function(hlay) {
        id <- unique(c(t(hlay))) 
        id[!is.na(id)]
    } 
    p <- grid.arrange(grobs = grob.list[select_grobs(lay)], layout_matrix = lay)

    print(p)
    return(NULL)
}

@transform_pandas(
    Output(rid="ri.vector.main.execute.8ac52849-bf13-40a9-b7f2-210f7fa33800"),
    corr_plot_d23vd22=Input(rid="ri.foundry.main.dataset.107ded67-aa74-4f7b-b01f-a950bc011336")
)
scatter_plot_v2_noline_skinny <- function(corr_plot_d23vd22) {
    #image: png
    imageWidth = 1300
    imageHeight = 1300
    dpi = 300

    png(
      filename=graphicsFile,
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")
          
    library(gridExtra)
    library(grid)
    library(ggplot2)
    library(lattice)
    library(stringr)
    library(ggpubr)

    df <- corr_plot_d23vd22 %>% select(-patient_id)

    colnames(df) <- gsub("-","", colnames(df))

    # IL6 (x) vs IP-10 (y)
    # TNFa (x) vs IP-10 (y)
    # TNFa (x) vs IL-6 (y)
    # IFNg (x) vs IP-10 (y)
    # IFNg (x) vs IL-6 (y)
    # IFNg (x) vs TNFa (y)
    # IL-15 (x) vs IP-10 (y)
    # IL-15 (x) vs IL-6 (y)
    # IL-15 (x) vs TNFa (y)
    # IL-15 (x) vs IFN (y)
    plots <- c("IL15-IFNg",
               "IL15-TNFa", "IFNg-TNFa",
               "IL15-IL6",  "IFNg-IL6",  "TNFa-IL6",
               "IL15-IP10", "IFNg-IP10", "TNFa-IP10", "IL6-IP10")
    lay <- rbind(c( 1, NA,  NA, NA),
                 c( 2,  3,  NA, NA),
                 c( 4,  5,   6, NA),
                 c( 7,  8,   9, 10))

    plot_cor <- function(name){
        text.size = 6

        x <- unlist(str_split(name,"-"))[1]
        y <- unlist(str_split(name,"-"))[2]

        if(x == "IL15"){
            axis.title.y = element_text(size = text.size+2, colour="black")
            axis.text.y=element_text(colour = "black", size=text.size)
        }else{
            axis.title.y = element_text(size = 2, colour="white")
            axis.text.y=element_text(colour = "white", size=1)
        }

        if(y == "IP10"){
            axis.title.x = element_text(size = text.size+2, colour="black")
            axis.text.x=element_text(colour="black", size=text.size)
        }else{
            axis.title.x = element_text(size = text.size+2, colour="white")
            axis.text.x=element_text(colour="white", size=text.size)
        }

        plt <- ggscatter(df, x=x, y=y, cor.method="spearman", color = "black", shape = 19, size = 0.2, # Points color, shape and size

                    conf.int = FALSE, # Add confidence interval
                    cor.coef = FALSE)+
                    theme(  axis.title.x = axis.title.x,
                            axis.title.y = axis.title.y,
                            axis.text.x = axis.text.x,
                            axis.text.y = axis.text.y,
                            axis.line = element_line(size=0.2),
                            axis.ticks= element_line(size=0.2))+
                    scale_y_continuous(labels = scales::number_format(accuracy = 0.1, trim=TRUE))+
                    scale_x_continuous(labels = scales::number_format(accuracy = 0.1, trim=TRUE))

        return(plt)
    }

    grob.list <- lapply(plots,plot_cor)
    select_grobs <- function(hlay) {
        id <- unique(c(t(hlay))) 
        id[!is.na(id)]
    } 
    p <- grid.arrange(grobs = grob.list[select_grobs(lay)], layout_matrix = lay)

    print(p)
    return(NULL)
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.42df9b3f-cd9d-4b74-81ef-96dd7f656f8b"),
    combined_complex_metadata=Input(rid="ri.foundry.main.dataset.8b5c5b80-01ad-4311-9156-3aad6c944a48"),
    merged_diff=Input(rid="ri.foundry.main.dataset.7ea602f6-f92c-4696-8f77-b5d39152978e")
)
vaccination_1I_filt_baseline <- function(merged_diff,combined_complex_metadata) {
    #This function uses pheatmap to draw a heatmap, scaling first by rows
    #(with samples in columns and genes in rows)
    # image: png
    suppressMessages(library(colorspace))
    suppressMessages(library(dendsort))
    suppressMessages(library(pheatmap))
    suppressMessages(library(tidyverse))
    suppressMessages(library(RColorBrewer))
    suppressMessages(library(dplyr))

    doheatmap <- function(dat, clus, clus2, ht, rn, cn, col) {
        require(pheatmap)
        require(dendsort)

        if (F) {
            tmean.scale = t(scale(t(dat)))
            tmean.scale = tmean.scale[!is.infinite(rowSums(tmean.scale)),]
            tmean.scale = na.omit(tmean.scale)
        } else {
            tmean.scale=dat
            #tmean.scale = na.omit(tmean.scale)
        }

        col_filter <- apply(as.data.frame(tmean.scale), 2, function(z) all(is.na(z)))
        row_filter <- apply(as.data.frame(tmean.scale), 1, function(z) length(z[is.na(z)]) > 10)
        tmean.scale <- tmean.scale[!row_filter, !col_filter]
        
        col.pal <- diverging_hcl(n=100, palette=col)
        if (FALSE) {
            col.pal = rev(col.pal)
        }
        # define metrics for clustering
        drows1 <- "correlation"
        dcols1 <- "correlation"

        minx = floor(min(tmean.scale))
        maxx = ceiling(max(tmean.scale))

        if (FALSE) {
            if(maxx > abs(minx)){
                minx = -1 * maxx
            }else{
                maxx = -1 * minx
            }
            breaks = seq(minx, maxx, length=100)
            legbreaks = seq(minx, maxx, length=5)
        } else {
            breaks = seq(-2.5, 2.5, length=100)
            legbreaks = seq(-2.5, 2.5, length=5)
        }
        breaks = sapply(breaks, signif, 4)
        legbreaks = sapply(legbreaks, signif, 4)

        if (TRUE) {
            treeheight <- 25
        } else {
            treeheight <- 0
        }

        gaps_col <- ""
        gaps_row <- ""

        gaps_col <- as.numeric(unlist(str_split(gaps_col,",")))
        gaps_row <- as.numeric(unlist(str_split(gaps_row,",")))

        gaps_col[is.na(gaps_col)] <- 0
        gaps_row[is.na(gaps_row)] <- 0

        hm.parameters <- list(
            tmean.scale, 
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
            gaps_col = gaps_col,
            gaps_row = gaps_row,
            annotation_names_col = FALSE,
            na_col = "#000000"
        )

        mat = t(tmean.scale)
        # print('calculated mat')
        
        do.call("pheatmap", c(hm.parameters))
    }

    df <- merged_diff

    samples_to_include = c("V002_d2-d1","V004_d2-d1","V005_d2-d1","V006_d2-d1","V007_d2-d1","V008_d2-d1","V009_d2-d1","V010_d2-d1","V011_d2-d1","V013_d2-d1","V015_d2-d1","V016_d2-d1","V017_d2-d1","V019_d2-d1","V020_d2-d1","V021_d2-d1","V022_d2-d1","V024_d2-d1","V025_d2-d1","V026_d2-d1","V027_d2-d1","V028_d2-d1","V029_d2-d1","V030_d2-d1","V031_d2-d1","V033_d2-d1","V034_d2-d1","V036_d2-d1","V037_d2-d1","V039_d2-d1","V041_d2-d1","V043_d2-d1","V045_d2-d1","V048_d2-d1","V049_d2-d1","V050_d2-d1","V051_d2-d1","V054_d2-d1","V056_d2-d1","V057_d2-d1","V058_d2-d1","V061_d2-d1","V063_d2-d1","V064_d2-d1","V066_d2-d1","V067_d2-d1","V068_d2-d1","V070_d2-d1","V071_d2-d1","V073_d2-d1","V077_d2-d1","V140_d2-d1","V141_d2-d1","V142_d2-d1","V143_d2-d1","V145_d2-d1","V147_d2-d1","V002_d23-d22","V004_d23-d22","V005_d23-d22","V006_d23-d22","V007_d23-d22","V008_d23-d22","V009_d23-d22","V011_d23-d22","V013_d23-d22","V015_d23-d22","V016_d23-d22","V017_d23-d22","V019_d23-d22","V020_d23-d22","V021_d23-d22","V022_d23-d22","V024_d23-d22","V025_d23-d22","V026_d23-d22","V027_d23-d22","V028_d23-d22","V029_d23-d22","V030_d23-d22","V031_d23-d22","V033_d23-d22","V034_d23-d22","V036_d23-d22","V037_d23-d22","V039_d23-d22","V041_d23-d22","V043_d23-d22","V045_d23-d22","V048_d23-d22","V049_d23-d22","V050_d23-d22","V051_d23-d22","V054_d23-d22","V056_d23-d22","V057_d23-d22","V058_d23-d22","V061_d23-d22","V063_d23-d22","V064_d23-d22","V066_d23-d22","V067_d23-d22","V068_d23-d22","V070_d23-d22","V071_d23-d22","V073_d23-d22","V077_d23-d22","V140_d23-d22","V141_d23-d22","V142_d23-d22","V143_d23-d22","V145_d23-d22","V147_d23-d22","V148_d23-d22","V002_d23-d22-d2-d1","V004_d23-d22-d2-d1","V005_d23-d22-d2-d1","V006_d23-d22-d2-d1","V007_d23-d22-d2-d1","V008_d23-d22-d2-d1","V009_d23-d22-d2-d1","V010_d23-d22-d2-d1","V011_d23-d22-d2-d1","V013_d23-d22-d2-d1","V015_d23-d22-d2-d1","V016_d23-d22-d2-d1","V017_d23-d22-d2-d1","V019_d23-d22-d2-d1","V020_d23-d22-d2-d1","V021_d23-d22-d2-d1","V022_d23-d22-d2-d1","V024_d23-d22-d2-d1","V025_d23-d22-d2-d1","V026_d23-d22-d2-d1","V027_d23-d22-d2-d1","V028_d23-d22-d2-d1","V029_d23-d22-d2-d1","V030_d23-d22-d2-d1","V031_d23-d22-d2-d1","V033_d23-d22-d2-d1","V034_d23-d22-d2-d1","V036_d23-d22-d2-d1","V037_d23-d22-d2-d1","V039_d23-d22-d2-d1","V041_d23-d22-d2-d1","V043_d23-d22-d2-d1","V045_d23-d22-d2-d1","V048_d23-d22-d2-d1","V049_d23-d22-d2-d1","V050_d23-d22-d2-d1","V051_d23-d22-d2-d1","V054_d23-d22-d2-d1","V056_d23-d22-d2-d1","V057_d23-d22-d2-d1","V058_d23-d22-d2-d1","V061_d23-d22-d2-d1","V063_d23-d22-d2-d1","V064_d23-d22-d2-d1")
    
    df.orig <- df %>% dplyr::select(one_of(c("Gene",samples_to_include)))
    df.mat = df.orig[ , (colnames(df.orig) != "Gene" )] %>% as.data.frame
    row.names(df.mat) <- df.orig$Gene
    df.mat <- as.data.frame(df.mat)
    
    annot <- combined_complex_metadata
    annot %>% dplyr::filter(sample_id %in% samples_to_include) -> annot
    annot$contrast <- factor(annot$contrast)
    print(annot$contrast)
    groups = c("contrast")
    relevel_factors <- FALSE
    factor_relevel <- c("contrast:d2_d1,d23_d22,vac2_vac1")
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
    #override colors
    colors <- c("#66C2A5","#984EA3","#B3B3B3")
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
      filename=graphicsFile,
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")

    p = doheatmap(dat=df.mat, clus=FALSE, clus2=TRUE, ht=50, rn=TRUE, cn=FALSE, col="Blue-Red 3")
    #Order output    
    row.order <- if("order" %in% names(p$tree_row)){
        p$tree_row[["order"]]
    }else{
        row.names(df.mat)
    }

    col.order <- if("order" %in% names(p$tree_col)){
        p$tree_col[["order"]]
    }else{
        colnames(df.mat)
    }

    df.mat <- df.mat[row.order,col.order]
    df.mat <- df.mat %>% rownames_to_column("Gene")

    return(df.mat)
}

@transform_pandas(
    Output(rid="ri.vector.main.execute.09c927dd-691d-46c8-a636-eba1d6ced62f"),
    combined_complex_metadata_immune=Input(rid="ri.foundry.main.dataset.7318826e-2f14-4395-8720-f2af0556bef7"),
    merged_diff_subset_immune=Input(rid="ri.foundry.main.dataset.a702eeb4-44cf-4dc4-85a5-462e19c6563f")
)
vaccination_1I_filt_baseline_immune <- function(merged_diff_subset_immune,combined_complex_metadata_immune) {
    #This function uses pheatmap to draw a heatmap, scaling first by rows
    #(with samples in columns and genes in rows)
    # image: png
    suppressMessages(library(colorspace))
    suppressMessages(library(dendsort))
    suppressMessages(library(pheatmap))
    suppressMessages(library(tidyverse))
    suppressMessages(library(RColorBrewer))
    suppressMessages(library(dplyr))

    doheatmap <- function(dat, clus, clus2, ht, rn, cn, col) {
        require(pheatmap)
        require(dendsort)

        if (F) {
            tmean.scale = t(scale(t(dat)))
            tmean.scale = tmean.scale[!is.infinite(rowSums(tmean.scale)),]
            tmean.scale = na.omit(tmean.scale)
        } else {
            tmean.scale=dat
            #tmean.scale = na.omit(tmean.scale)
        }

        col_filter <- apply(as.data.frame(tmean.scale), 2, function(z) all(is.na(z)))
        row_filter <- apply(as.data.frame(tmean.scale), 1, function(z) length(z[is.na(z)]) > 5)
        tmean.scale <- tmean.scale[!row_filter, !col_filter]

        col.pal <- diverging_hcl(n=100, palette=col)
        if (FALSE) {
            col.pal = rev(col.pal)
        }
        # define metrics for clustering
        drows1 <- "correlation"
        dcols1 <- "correlation"

        minx = floor(min(tmean.scale))
        maxx = ceiling(max(tmean.scale))

        if (FALSE) {
            if(maxx > abs(minx)){
                minx = -1 * maxx
            }else{
                maxx = -1 * minx
            }
            breaks = seq(minx, maxx, length=100)
            legbreaks = seq(minx, maxx, length=5)
        } else {
            breaks = seq(-3.5, 3.5, length=100)
            legbreaks = seq(-3.5, 3.5, length=5)
        }
        breaks = sapply(breaks, signif, 4)
        legbreaks = sapply(legbreaks, signif, 4)

        if (TRUE) {
            treeheight <- 25
        } else {
            treeheight <- 0
        }

        gaps_col <- ""
        gaps_row <- ""

        gaps_col <- as.numeric(unlist(str_split(gaps_col,",")))
        gaps_row <- as.numeric(unlist(str_split(gaps_row,",")))

        gaps_col[is.na(gaps_col)] <- 0
        gaps_row[is.na(gaps_row)] <- 0

        hm.parameters <- list(
            tmean.scale, 
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
            gaps_col = gaps_col,
            gaps_row = gaps_row,
            annotation_names_col = FALSE,
            na_col = "#000000"
        )

        mat = t(tmean.scale)
        # print('calculated mat')
        
        do.call("pheatmap", c(hm.parameters))
    }

    df <- merged_diff_subset_immune

    samples_to_include = c("V001_d2-d1","V003_d2-d1","V014_d2-d1","V032_d2-d1","V062_d2-d1","V001_d23-d22","V003_d23-d22","V014_d23-d22","V032_d23-d22","V062_d23-d22","V001_d23-d22-d2-d1","V003_d23-d22-d2-d1","V014_d23-d22-d2-d1","V032_d23-d22-d2-d1","V062_d23-d22-d2-d1")

    df.orig <- df %>% dplyr::select(one_of(c("Gene",samples_to_include)))
    df.mat = df.orig[ , (colnames(df.orig) != "Gene" )] %>% as.data.frame
    row.names(df.mat) <- df.orig$Gene
    df.mat <- as.data.frame(df.mat)
    
    annot <- combined_complex_metadata_immune
    annot %>% dplyr::filter(sample_id %in% samples_to_include) -> annot
    annot$contrast <- factor(annot$contrast)

    groups = c("contrast")
    relevel_factors <- FALSE
    factor_relevel <- c("contrast:d2_d1,d23_d22,vac2_vac1")
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
    colors <- c("#66C2A5","#984EA3","#B3B3B3")

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
      filename=graphicsFile,
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")

    p = doheatmap(dat=df.mat, clus=FALSE, clus2=TRUE, ht=50, rn=TRUE, cn=FALSE, col="Blue-Red 3")
    #Order output    
    row.order <- if("order" %in% names(p$tree_row)){
        p$tree_row[["order"]]
    }else{
        row.names(df.mat)
    }

    col.order <- if("order" %in% names(p$tree_col)){
        p$tree_col[["order"]]
    }else{
        colnames(df.mat)
    }

    df.mat <- df.mat[row.order,col.order]
    df.mat <- df.mat %>% rownames_to_column("Gene")

    return(df.mat)
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.9e0f4dae-b081-4814-836f-585369a7c3aa"),
    naive_deg_analysis=Input(rid="ri.foundry.main.dataset.862f8327-dcb8-414e-b8ae-fb74aaa3a15b")
)
volcano_IIvI <- function(naive_deg_analysis) {
    # image: png

    # Changelog
    # 2020-10-29 Add support for pval == 0

    suppressMessages(library(stringr))
    suppressMessages(library(tidyverse))
    suppressMessages(library(ggplot2))
    
    #Class handling
    df <- naive_deg_analysis
    if(class(df) == "SparkDataFrame"){
        df <- SparkR::collect(df)
    }
  
    # Log fold-change column
    # p-value column
    # Position of p-value threshold line #default 0.05
    # Position of log fold-change threshold lines #default 1 lfc
    # Label column #default gene names
    # Specific genes to label
    # Number of genes to label
    
    #label topN genes or specific genes
    
    label.col <- "Gene"
    sig.col <- "d23-d22-d2-d1_adjpval"
    lfc.col <- "d23-d22-d2-d1_logFC"
    columns_of_interest <- c(label.col,lfc.col,sig.col)

    df <- df %>% dplyr::select(one_of(columns_of_interest)) %>% dplyr::filter(!is.na(!!rlang::sym(lfc.col))) #seurat introduces NAs

    change_lfc_name <- "log2FC"
    change_sig_name <- "FDR"
    colnames(df) <- c(label.col,change_lfc_name,change_sig_name)

    cat(paste0("Genes in initial dataset: ", nrow(df),"\n"))

    #Select top genes by logFC or Significance
    no_genes_to_label <- 30
    value_to_sort_the_output_dataset <- "p-value"
    if (value_to_sort_the_output_dataset=="fold-change") {
        df %>% dplyr::arrange(desc(abs(!!rlang::sym(change_lfc_name)))) -> df
    } else if (value_to_sort_the_output_dataset=="p-value") {
        df %>% dplyr::arrange(!!rlang::sym(change_sig_name)) -> df
    }
    genes_to_label <- as.character(df[1:no_genes_to_label,label.col])

    additional_labels <- ""
    additional_labels <- unlist(str_split(additional_labels,","))
    filter <- additional_labels %in% df[,label.col]
    additional_labels <- additional_labels[filter]
    missing_labels <- additional_labels[!filter]

    if(!is.na(missing_labels)){
        cat("Could not find:\n")
        print(missing_labels)
    }
    use_only_addition_labels <- FALSE
    if(use_only_addition_labels){
        genes_to_label <- additional_labels
    }else{
        genes_to_label <- unique(append(genes_to_label, additional_labels))
    }

    # Plot variables
    title <- "Volcano Plots"
    subtitle <- "Vaccination 2 vs Vaccination 1 (1 dpv vs pre)"

    pCutoff  = 0.05
    FCcutoff = 1.0

    significant=as.vector(table(abs( df[,change_lfc_name] ) > FCcutoff &
                                     df[,change_sig_name]   < pCutoff))[2]

    # fix pvalue == 0
    shapeCustom <- rep(19,nrow(df))
    maxy <-  max(-log10(df[[change_sig_name]]), na.rm=TRUE)
    if(0 > 0){
        maxy <- 0
    }
    
    cat(paste0("Maxy: ",maxy,"\n"))
    if(maxy == Inf){
        # Sometimes, pvalues == 0
        keep <- df[[change_sig_name]] > 0
        df[[change_sig_name]][!keep] <- min(df[[change_sig_name]][keep])
        shapeCustom[!keep] <- 17

        maxy <- -log10(min(df[[change_sig_name]][keep]))
        cat("Some p-values equal zero. Adjusting y-limits.\n")
        cat(paste0("Maxy adjusted: ",maxy,"\n"))

    }

    # By default, nothing will be greater than maxy. User can set this value lower
    keep <- -log10(df[[change_sig_name]]) <= maxy
    df[[change_sig_name]][!keep] <- maxy
    shapeCustom[!keep] <- 17

    names(shapeCustom)<- rep("Exact",length(shapeCustom))
    names(shapeCustom)[shapeCustom == 17] <- "Adjusted"
    
    #Remove if nothin' doin'
    if(all(shapeCustom == 19)){
        shapeCustom <- NULL
    }
    
    maxy <- ceiling(maxy)

    use_custom_xlab <- FALSE
    if(use_custom_xlab){
        xlab <- gsub("_"," ",change_lfc_name)
    }else{
        xlab <- bquote(~Log[2]~ "fold change")
    }

    xlim_additional <- 0
    ylim_additional <- 0
    axisLabSize <- 24
    labSize <- 4
    pointSize <- 2

    imageWidth = 3000
    imageHeight = 3000
    dpi = 300

    png(
      filename=graphicsFile,
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")

    p <- EnhancedVolcano(   df,x=change_lfc_name,y=change_sig_name,
                            lab=df[,label.col],
                            selectLab = genes_to_label,
                            title=paste0(title," (Significant=",significant,")"),
                            subtitle <- subtitle,
                            xlab=xlab,
                            ylab=bquote(~-Log[10]~.(change_sig_name)),
                            xlim=c(floor(min(df[,change_lfc_name])) - xlim_additional,ceiling(max(df[,change_lfc_name]))+ xlim_additional),
                            ylim=c(0, maxy + ylim_additional),
                            pCutoff=pCutoff,
                            FCcutoff=FCcutoff,
                            axisLabSize=axisLabSize,
                            labSize=labSize,
                            pointSize=pointSize,
                            shapeCustom=shapeCustom
                            )
    print(p)

    df$rank <- -log10(df[[change_sig_name]]) * sign(df[[change_lfc_name]]) #sig already -log10
    df <- df %>% arrange(desc(rank))
    return(df)
}

# This is directly copied from https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html v1.6.0
EnhancedVolcano <- function(
  toptable,
  lab,
  x,
  y,
  selectLab = NULL,
  xlim = c(min(toptable[[x]], na.rm=TRUE) - 1.5,
    max(toptable[[x]], na.rm=TRUE) + 1.5),
  ylim = c(0, max(-log10(toptable[[y]]), na.rm=TRUE) + 5),
  xlab = bquote(~Log[2]~ "fold change"),
  ylab = bquote(~-Log[10]~italic(P)),
  axisLabSize = 18,
  title = 'Volcano plot',
  subtitle = bquote(italic(EnhancedVolcano)),
  caption = paste0('total = ', nrow(toptable), ' variables'),
  titleLabSize = 18,
  subtitleLabSize = 14,
  captionLabSize = 14,
  pCutoff = 10e-6,
  FCcutoff = 1.0,
  cutoffLineType = 'longdash',
  cutoffLineCol = 'black',
  cutoffLineWidth = 0.4,
  pointSize = 2.0,
  labSize = 5.0,
  labCol = 'black',
  labFace = 'plain',
  labhjust = 0.5,
  labvjust = 1.5,
  boxedLabels = FALSE,
  shape = 19,
  shapeCustom = NULL,
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  colCustom = NULL,
  colAlpha = 1/2,
  colGradient = NULL,
  colGradientBreaks = c(pCutoff, 1.0),
  colGradientLabels = c('0', '1.0'),
  colGradientLimits = c(0, 1.0),
  legendLabels = c('NS', expression(Log[2]~FC),
    'p-value', expression(p-value~and~log[2]~FC)),
  legendPosition = 'top',
  legendLabSize = 14,
  legendIconSize = 5.0,
  legendDropLevels = TRUE,
  encircle = NULL,
  encircleCol = 'black',
  encircleFill = 'pink',
  encircleAlpha = 3/4,
  encircleSize = 2.5,
  shade = NULL,
  shadeFill = 'grey',
  shadeAlpha = 1/2,
  shadeSize = 0.01,
  shadeBins = 2,
  drawConnectors = FALSE,
  widthConnectors = 0.5,
  typeConnectors = 'closed',
  endsConnectors = 'first',
  lengthConnectors = unit(0.01, 'npc'),
  colConnectors = 'grey10',
  arrowheads = TRUE,
  hline = NULL,
  hlineType = 'longdash',
  hlineCol = 'black',
  hlineWidth = 0.4,
  vline = NULL,
  vlineType = 'longdash',
  vlineCol = 'black',
  vlineWidth = 0.4,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  border = 'partial',
  borderWidth = 0.8,
  borderColour = 'black', 
  raster = FALSE)
{
  if(!is.numeric(toptable[[x]])) {
    stop(paste(x, ' is not numeric!', sep=''))
  }

  if(!is.numeric(toptable[[y]])) {
    stop(paste(y, ' is not numeric!', sep=''))
  }
  
  if (raster) {
    geom_point <- geom_point_rast
  }

  i <- xvals <- yvals <- Sig <- NULL

  toptable <- as.data.frame(toptable)
  toptable$Sig <- 'NS'
  toptable$Sig[(abs(toptable[[x]]) > FCcutoff)] <- 'FC'
  toptable$Sig[(toptable[[y]] < pCutoff)] <- 'P'
  toptable$Sig[(toptable[[y]] < pCutoff) &
    (abs(toptable[[x]]) > FCcutoff)] <- 'FC_P'
  toptable$Sig <- factor(toptable$Sig,
    levels=c('NS','FC','P','FC_P'))

  # some software programs return 0 for very low p-values
  # These throw an error in EnhancedVolcano
  # Detect these, issue warning, and convert these to
  # machine-lowest value possible
  #####
  # New functionality in > v1.2:
  # Now convert to 10^-1 lower than lowest non-zero p-value
  if (min(toptable[[y]], na.rm=TRUE) == 0) {
    # <= v1.2
    #warning(paste("One or more P values is 0.",
    #  "Converting to minimum possible value..."),
    #  call. = FALSE)
    #toptable[which(toptable[[y]] == 0), y] <- .Machine$double.xmin
    warning(paste('One or more p-values is 0.',
      'Converting to 10^-1 * current',
      'lowest non-zero p-value...'),
      call. = FALSE)
    toptable[which(toptable[[y]] == 0), y] <- min(
      toptable[which(toptable[[y]] != 0), y],
      na.rm = TRUE) * 10^-1
  }

  toptable$lab <- lab
  toptable$xvals <- toptable[[x]]
  toptable$yvals <- toptable[[y]]

  # If user has supplied values in selectLab, convert labels to
  # NA and then re-set with those in selectLab
  if (!is.null(selectLab)) {
    names.new <- rep(NA, length(toptable$lab))
    indices <- which(toptable$lab %in% selectLab)
    names.new[indices] <- toptable$lab[indices]
    toptable$lab <- names.new
  }

  # create a base theme that will later be modified
  th <- theme_bw(base_size = 24) +

    theme(
      legend.background = element_rect(),

      # title, subtitle, and caption
      plot.title = element_text(
        angle = 0,
        size = titleLabSize,
        face = 'bold',
        vjust = 1),
      plot.subtitle = element_text(
        angle = 0,
        size = subtitleLabSize,
        face = 'plain',
        vjust = 1),
      plot.caption = element_text(
        angle = 0,
        size = captionLabSize,
        face = 'plain',
        vjust = 1),

      # axis text
      axis.text.x = element_text(
        angle = 0,
        size = axisLabSize,
        vjust = 1),
      axis.text.y = element_text(
        angle = 0,
        size = axisLabSize,
        vjust = 0.5),
      axis.title = element_text(
        size = axisLabSize),

      # legend
      legend.position = legendPosition,
      legend.key = element_blank(),
      legend.key.size = unit(0.5, 'cm'),
      legend.text = element_text(
        size = legendLabSize),
      title = element_text(
        size = legendLabSize),
      legend.title = element_blank())

  # Create the plot object differently based on whether colCustom 
  # and shapeCustom are NULL or not. This helps to avoid messing up
  # the legend.
  #
  # 1, both colCustom and shapeCustom are activated
  if (!is.null(colCustom) & !is.null(shapeCustom)) {

    plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

      # over-ride legend icon sizes for colour and shape.
      # guide_legends are separate for colour and shape;
      # so, legends will be drawn separate
      guides(
        colour = guide_legend(
          order = 1,
          override.aes = list(
            size = legendIconSize)),
        shape = guide_legend(
          order = 2,
          override.aes = list(
            size = legendIconSize))) +

      # include new shape and colour encodings as aes
      geom_point(
        aes(
          color = factor(names(colCustom)),
          shape = factor(names(shapeCustom))),
        alpha = colAlpha,
        size = pointSize,
        na.rm = TRUE) +

      # specify the colour and shape with the supplied encoding
      scale_color_manual(values = colCustom) +
      scale_shape_manual(values = shapeCustom)

  # 2, only colCustom is activated and 'shape' has just a single value
  } else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 1) {

    plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

      # over-ride legend icon sizes for colour and shape.
      # guide_legends are separate for colour and shape;
      # so, legends will be drawn separate IF shape is also
      # included as aes to geom_point (it is not, here)
      guides(
        colour = guide_legend(
          order = 1,
          override.aes = list(
            size = legendIconSize)),
        shape = guide_legend(
          order = 2,
          override.aes = list(
            size = legendIconSize))) +

      # include new colour encodings as aes.
      # 'shape' is included, but outside aes
      geom_point(
        aes(
          color = factor(names(colCustom))),
        alpha = colAlpha,
        shape = shape,
        size = pointSize,
        na.rm = TRUE) +

      # specify the colour with the supplied encoding
      scale_color_manual(values = colCustom) +

      # 'shape' is not included as aes. Specifying guide = TRUE
      # here will result in legends merging
      scale_shape_manual(guide = TRUE)

  # 3, only colCustom is activated and 'shape' has 4 values
  } else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 4) {

    plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

      # over-ride legend icon sizes for colour and shape.
      # guide_legends are separate for colour and shape;
      # so, legends will be drawn separate
      guides(
        colour = guide_legend(
          order = 1,
          override.aes = list(
            size = legendIconSize)),
        shape = guide_legend(
          order = 2,
          override.aes = list(
            size = legendIconSize))) +

      # include new colour encodings as aes.
      # 'shape' is included in aes and mapped to 4
      # categories of NS, FC, P, FC_P
      geom_point(
        aes(
          color = factor(names(colCustom)),
          shape = Sig),
        alpha = colAlpha,
        size = pointSize,
        na.rm = TRUE) +

      # specify the colour with the supplied encoding
      scale_color_manual(values = colCustom) +

      # as it is included as aes, a separate legend
      # for 'shape' will be drawn. Here, over-ride that
      # legend
      scale_shape_manual(
        values = c(
          NS = shape[1],
          FC = shape[2],
          P = shape[3],
          FC_P = shape[4]),
        labels = c(
          NS = legendLabels[1],
          FC = legendLabels[2],
          P = legendLabels[3],
          FC_P = legendLabels[4]),
        guide = TRUE,
        drop = legendDropLevels)

  # 4, only shapeCustom is activated
  } else if (is.null(colCustom) & !is.null(shapeCustom)) {

    if (is.null(colGradient)) {

      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        # over-ride legend icon sizes for colour and shape.
        # guide_legends are separate for colour and shape;
        # so, legends will be drawn separate
        guides(
          colour = guide_legend(
            order = 1,
            override.aes = list(
              size = legendIconSize)),
          shape = guide_legend(
            order = 2,
            override.aes = list(
              size = legendIconSize))) +

        # include new shape encodings as aes.
        # Standard colour for NS, FC, P, FC_P,
        # are added to aes, too.
        geom_point(
          aes(
            color = Sig,
            shape = factor(names(shapeCustom))),
          alpha = colAlpha,
          size = pointSize,
          na.rm = TRUE) +

        # as it is included as aes, a separate legend
        # for 'colour' will be drawn. Here, over-ride that
        # legend
        scale_color_manual(
          values = c(
            NS = col[1],
            FC = col[2],
            P = col[3],
            FC_P = col[4]),
          labels = c(
            NS = legendLabels[1],
            FC = legendLabels[2],
            P = legendLabels[3],
            FC_P = legendLabels[4]),
          drop = legendDropLevels) +

        # specify the shape with the supplied encoding
        scale_shape_manual(values = shapeCustom)

    } else {

      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        # over-ride legend icon sizes for colour and shape.
        # guide_legends are separate for colour and shape;
        # so, legends will be drawn separate
        guides(
          shape = guide_legend(
            order = 2,
            override.aes = list(
              size = legendIconSize))) +

        # include new shape encodings as aes.
        # Standard colour for NS, FC, P, FC_P,
        # are added to aes, too.
        geom_point(
          aes(
            color = Sig,
            shape = factor(names(shapeCustom))),
          alpha = colAlpha,
          size = pointSize,
          na.rm = TRUE) +

        scale_colour_gradient(
          low = colGradient[1],
          high = colGradient[2],
          limits = colGradientLimits,
          breaks = colGradientBreaks,
          labels = colGradientLabels)

        # specify the shape with the supplied encoding
        scale_shape_manual(values = shapeCustom)

    }

  # 5, both colCustom and shapeCustom are null;
  # only a single shape value specified
  } else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 1) {

    if (is.null(colGradient)) {

      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        # over-ride legend icon sizes for colour and shape.
        # including 'shape' in the colour guide_legend here
        # results in the legends merging
        guides(colour = guide_legend(
          order = 1,
          override.aes = list(
            shape = shape,
            size = legendIconSize))) +

        geom_point(
          aes(color = Sig),
          alpha = colAlpha,
          shape = shape,
          size = pointSize,
          na.rm = TRUE) +

        scale_color_manual(
          values = c(
            NS = col[1],
            FC = col[2],
            P = col[3],
            FC_P = col[4]),
          labels = c(
            NS = legendLabels[1],
            FC = legendLabels[2],
            P = legendLabels[3],
            FC_P = legendLabels[4]),
          drop = legendDropLevels)

    } else {

      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        geom_point(
          aes(color = yvals),
          alpha = colAlpha,
          shape = shape,
          size = pointSize,
          na.rm = TRUE) +

        scale_colour_gradient(
          low = colGradient[1],
          high = colGradient[2],
          limits = colGradientLimits,
          breaks = colGradientBreaks,
          labels = colGradientLabels)
    }

  # 6, both colCustom and shapeCustom are null;
  # four shape values are specified
  } else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 4) {

    if (is.null(colGradient)) {
      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        # over-ride legend icon sizes for colour and shape.
        # including 'shape' in the colour guide_legend here
        # results in the legends merging
        guides(colour = guide_legend(
          order = 1,
          override.aes = list(
            shape = c(
              NS = shape[1],
              FC = shape[2],
              P = shape[3],
              FC_P = shape[4]),
            size = legendIconSize))) +

        geom_point(
          aes(
            color = Sig,
            shape = Sig),
          alpha = colAlpha,
          size = pointSize,
          na.rm = TRUE) +

        scale_color_manual(
          values = c(
            NS = col[1],
            FC = col[2],
            P = col[3],
            FC_P = col[4]),
          labels = c(
            NS = legendLabels[1],
            FC = legendLabels[2],
            P = legendLabels[3],
            FC_P = legendLabels[4]),
          drop = legendDropLevels) +

        scale_shape_manual(
          values = c(
            NS = shape[1],
            FC = shape[2],
            P = shape[3],
            FC_P = shape[4]),
          guide = FALSE,
          drop = legendDropLevels)

    } else {

      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        geom_point(
          aes(
            color = yvals,
            shape = Sig),
          alpha = colAlpha,
          size = pointSize,
          na.rm = TRUE) +

        scale_colour_gradient(
          low = colGradient[1],
          high = colGradient[2],
          limits = colGradientLimits,
          breaks = colGradientBreaks,
          labels = colGradientLabels) +

        scale_shape_manual(
          values = c(
            NS = shape[1],
            FC = shape[2],
            P = shape[3],
            FC_P = shape[4]),
          guide = FALSE,
          drop = legendDropLevels)

    }
  }

  # add more elements to the plot
  plot <- plot +

    xlab(xlab) +
    ylab(ylab) +

    xlim(xlim[1], xlim[2]) +
    ylim(ylim[1], ylim[2]) +

    geom_vline(xintercept = c(-FCcutoff, FCcutoff),
      linetype = cutoffLineType,
      colour = cutoffLineCol,
      size = cutoffLineWidth) +

    geom_hline(yintercept = -log10(pCutoff),
      linetype = cutoffLineType,
      colour = cutoffLineCol,
      size = cutoffLineWidth)

  # add elements to the plot for title, subtitle, caption
  plot <- plot + labs(title = title, 
    subtitle = subtitle, caption = caption)

  # add elements to the plot for vlines and hlines
  if (!is.null(vline)) {
    plot <- plot + geom_vline(xintercept = vline,
      linetype = vlineType,
      colour = vlineCol,
      size = vlineWidth)
  }
  if (!is.null(hline)) {
    plot <- plot + geom_hline(yintercept = -log10(hline),
      linetype = hlineType,
      colour = hlineCol,
      size = hlineWidth)
  }

  # Border around plot
  if (border == 'full') {
    plot <- plot + theme(panel.border = element_rect(
      colour = borderColour, fill = NA, size = borderWidth))
  } else if (border == 'partial') {
    plot <- plot + theme(axis.line = element_line(
      size = borderWidth, colour = borderColour),
      panel.border = element_blank(),
      panel.background = element_blank())
  } else {
    stop('Unrecognised value passed to \'border\'. Must be \'full\' or \'partial\'')
  }

  # Gridlines
  if (gridlines.major) {
    plot <- plot + theme(panel.grid.major = element_line())
  } else {
    plot <- plot + theme(panel.grid.major = element_blank())
  }
  if (gridlines.minor) {
    plot <- plot + theme(panel.grid.minor = element_line())
  } else {
    plot <- plot + theme(panel.grid.minor = element_blank())
  }

  # user has specified to draw with geom_text or geom_label?
  if (!boxedLabels) {

    # For labeling with geom_[text|label]_repel and
    # geom_[text|label] with check_overlap = TRUE, 4 possible
    # scenarios can arise
    if (drawConnectors && is.null(selectLab)) {

      if (arrowheads) {
        arr <- arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors)
      } else {
        arr <- NULL
      }

      plot <- plot + geom_text_repel(
        data = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff),
        aes(label = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff)[["lab"]]),
        size = labSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arr,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    } else if (drawConnectors && !is.null(selectLab)) {

      if (arrowheads) {
        arr <- arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors)
      } else {
        arr <- NULL
      }

      plot <- plot + geom_text_repel(
        data = subset(toptable,
          !is.na(toptable[['lab']])),
        aes(label = subset(toptable,
          !is.na(toptable[['lab']]))[['lab']]),
        size = labSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arr,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    } else if (!drawConnectors && !is.null(selectLab)) {

      plot <- plot + geom_text(
        data = subset(toptable,
          !is.na(toptable[['lab']])),
        aes(
          label = subset(toptable,
            !is.na(toptable[['lab']]))[['lab']]),
        size = labSize,
        check_overlap = TRUE,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    } else if (!drawConnectors && is.null(selectLab)) {

      plot <- plot + geom_text(
        data = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff),
        aes(label = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff)[['lab']]),
        size = labSize,
        check_overlap = TRUE,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)
    }

  } else {

    # For labeling with geom_[text|label]_repel and
    # geom_[text|label] with check_overlap = TRUE, 4 possible
    # scenarios can arise
    if (drawConnectors && is.null(selectLab)) {

      if (arrowheads) {
        arr <- arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors)
      } else {
        arr <- NULL
      }

      plot <- plot + geom_label_repel(
        data = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff),
        aes(label = subset(toptable,
          toptable[[y]]<pCutoff &
            abs(toptable[[x]]) > FCcutoff)[['lab']]),
        size = labSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arr,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    } else if (drawConnectors && !is.null(selectLab)) {

      if (arrowheads) {
        arr <- arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors)
      } else {
        arr <- NULL
      }

      plot <- plot + geom_label_repel(
        data = subset(toptable,
          !is.na(toptable[['lab']])),
        aes(label = subset(toptable,
          !is.na(toptable[['lab']]))[['lab']]),
        size = labSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arr,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    } else if (!drawConnectors && !is.null(selectLab)) {

      plot <- plot + geom_label(
        data = subset(toptable,
          !is.na(toptable[["lab"]])),
        aes(
          label = subset(toptable,
            !is.na(toptable[['lab']]))[['lab']]),
        size = labSize,
        #check_overlap = TRUE,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    } else if (!drawConnectors && is.null(selectLab)) {

      plot <- plot + geom_label(
        data = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff),
        aes(label = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff)[['lab']]),
        size = labSize,
        #check_overlap = TRUE,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    }
  }

  # encircle
  if (!is.null(encircle)) {
    plot <- plot + 
      geom_encircle(
        data = subset(toptable,
          rownames(toptable) %in% encircle),
        colour = encircleCol,
        fill = encircleFill,
        alpha = encircleAlpha,
        size = encircleSize,
        show.legend = FALSE,
        na.rm = TRUE)
  }

  # shade
  if (!is.null(shade)) {
    plot <- plot + 
      stat_density2d(
        data = subset(toptable,
          rownames(toptable) %in% shade),
        fill = shadeFill,
        alpha = shadeAlpha,
        geom = 'polygon',
        contour = TRUE,
        size = shadeSize,
        bins = shadeBins,
        show.legend = FALSE,
        na.rm = TRUE)
  }

  return(plot)
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.b2f0eb24-a391-4946-a354-963e757c480b"),
    naive_deg_analysis=Input(rid="ri.foundry.main.dataset.862f8327-dcb8-414e-b8ae-fb74aaa3a15b")
)
volcano_d23vd22 <- function(naive_deg_analysis) {
    # image: png

    # Changelog
    # 2020-10-29 Add support for pval == 0

    suppressMessages(library(stringr))
    suppressMessages(library(tidyverse))
    suppressMessages(library(ggplot2))
    
    #Class handling
    df <- naive_deg_analysis
    if(class(df) == "SparkDataFrame"){
        df <- SparkR::collect(df)
    }
  
    # Log fold-change column
    # p-value column
    # Position of p-value threshold line #default 0.05
    # Position of log fold-change threshold lines #default 1 lfc
    # Label column #default gene names
    # Specific genes to label
    # Number of genes to label
    
    #label topN genes or specific genes
    
    label.col <- "Gene"
    sig.col <- "d23-d22_adjpval"
    lfc.col <- "d23-d22_logFC"
    columns_of_interest <- c(label.col,lfc.col,sig.col)

    df <- df %>% dplyr::select(one_of(columns_of_interest)) %>% dplyr::filter(!is.na(!!rlang::sym(lfc.col))) #seurat introduces NAs

    change_lfc_name <- "log2FC"
    change_sig_name <- "FDR"
    colnames(df) <- c(label.col,change_lfc_name,change_sig_name)

    cat(paste0("Genes in initial dataset: ", nrow(df),"\n"))

    #Select top genes by logFC or Significance
    no_genes_to_label <- 30
    value_to_sort_the_output_dataset <- "p-value"
    if (value_to_sort_the_output_dataset=="fold-change") {
        df %>% dplyr::arrange(desc(abs(!!rlang::sym(change_lfc_name)))) -> df
    } else if (value_to_sort_the_output_dataset=="p-value") {
        df %>% dplyr::arrange(!!rlang::sym(change_sig_name)) -> df
    }
    genes_to_label <- as.character(df[1:no_genes_to_label,label.col])

    additional_labels <- ""
    additional_labels <- unlist(str_split(additional_labels,","))
    filter <- additional_labels %in% df[,label.col]
    additional_labels <- additional_labels[filter]
    missing_labels <- additional_labels[!filter]

    if(!is.na(missing_labels)){
        cat("Could not find:\n")
        print(missing_labels)
    }
    use_only_addition_labels <- FALSE
    if(use_only_addition_labels){
        genes_to_label <- additional_labels
    }else{
        genes_to_label <- unique(append(genes_to_label, additional_labels))
    }

    # Plot variables
    title <- "Volcano Plots"
    subtitle <- "Vaccination 2 ( 1 dpv vs pre)"

    pCutoff  = 0.05
    FCcutoff = 1.0

    significant=as.vector(table(abs( df[,change_lfc_name] ) > FCcutoff &
                                     df[,change_sig_name]   < pCutoff))[2]

    # fix pvalue == 0
    shapeCustom <- rep(19,nrow(df))
    maxy <-  max(-log10(df[[change_sig_name]]), na.rm=TRUE)
    if(0 > 0){
        maxy <- 0
    }
    
    cat(paste0("Maxy: ",maxy,"\n"))
    if(maxy == Inf){
        # Sometimes, pvalues == 0
        keep <- df[[change_sig_name]] > 0
        df[[change_sig_name]][!keep] <- min(df[[change_sig_name]][keep])
        shapeCustom[!keep] <- 17

        maxy <- -log10(min(df[[change_sig_name]][keep]))
        cat("Some p-values equal zero. Adjusting y-limits.\n")
        cat(paste0("Maxy adjusted: ",maxy,"\n"))

    }

    # By default, nothing will be greater than maxy. User can set this value lower
    keep <- -log10(df[[change_sig_name]]) <= maxy
    df[[change_sig_name]][!keep] <- maxy
    shapeCustom[!keep] <- 17

    names(shapeCustom)<- rep("Exact",length(shapeCustom))
    names(shapeCustom)[shapeCustom == 17] <- "Adjusted"
    
    #Remove if nothin' doin'
    if(all(shapeCustom == 19)){
        shapeCustom <- NULL
    }
    
    maxy <- ceiling(maxy)

    use_custom_xlab <- FALSE
    if(use_custom_xlab){
        xlab <- gsub("_"," ",change_lfc_name)
    }else{
        xlab <- bquote(~Log[2]~ "fold change")
    }

    xlim_additional <- 0
    ylim_additional <- 0
    axisLabSize <- 24
    labSize <- 4
    pointSize <- 2

    imageWidth = 3000
    imageHeight = 3000
    dpi = 300

    png(
      filename=graphicsFile,
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")

    p <- EnhancedVolcano(   df,x=change_lfc_name,y=change_sig_name,
                            lab=df[,label.col],
                            selectLab = genes_to_label,
                            title=paste0(title," (Significant=",significant,")"),
                            subtitle <- subtitle,
                            xlab=xlab,
                            ylab=bquote(~-Log[10]~.(change_sig_name)),
                            xlim=c(floor(min(df[,change_lfc_name])) - xlim_additional,ceiling(max(df[,change_lfc_name]))+ xlim_additional),
                            ylim=c(0, maxy + ylim_additional),
                            pCutoff=pCutoff,
                            FCcutoff=FCcutoff,
                            axisLabSize=axisLabSize,
                            labSize=labSize,
                            pointSize=pointSize,
                            shapeCustom=shapeCustom
                            )
    print(p)

    df$rank <- -log10(df[[change_sig_name]]) * sign(df[[change_lfc_name]]) #sig already -log10
    df <- df %>% arrange(desc(rank))
    return(df)
}

# This is directly copied from https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html v1.6.0
EnhancedVolcano <- function(
  toptable,
  lab,
  x,
  y,
  selectLab = NULL,
  xlim = c(min(toptable[[x]], na.rm=TRUE) - 1.5,
    max(toptable[[x]], na.rm=TRUE) + 1.5),
  ylim = c(0, max(-log10(toptable[[y]]), na.rm=TRUE) + 5),
  xlab = bquote(~Log[2]~ "fold change"),
  ylab = bquote(~-Log[10]~italic(P)),
  axisLabSize = 18,
  title = 'Volcano plot',
  subtitle = bquote(italic(EnhancedVolcano)),
  caption = paste0('total = ', nrow(toptable), ' variables'),
  titleLabSize = 18,
  subtitleLabSize = 14,
  captionLabSize = 14,
  pCutoff = 10e-6,
  FCcutoff = 1.0,
  cutoffLineType = 'longdash',
  cutoffLineCol = 'black',
  cutoffLineWidth = 0.4,
  pointSize = 2.0,
  labSize = 5.0,
  labCol = 'black',
  labFace = 'plain',
  labhjust = 0.5,
  labvjust = 1.5,
  boxedLabels = FALSE,
  shape = 19,
  shapeCustom = NULL,
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  colCustom = NULL,
  colAlpha = 1/2,
  colGradient = NULL,
  colGradientBreaks = c(pCutoff, 1.0),
  colGradientLabels = c('0', '1.0'),
  colGradientLimits = c(0, 1.0),
  legendLabels = c('NS', expression(Log[2]~FC),
    'p-value', expression(p-value~and~log[2]~FC)),
  legendPosition = 'top',
  legendLabSize = 14,
  legendIconSize = 5.0,
  legendDropLevels = TRUE,
  encircle = NULL,
  encircleCol = 'black',
  encircleFill = 'pink',
  encircleAlpha = 3/4,
  encircleSize = 2.5,
  shade = NULL,
  shadeFill = 'grey',
  shadeAlpha = 1/2,
  shadeSize = 0.01,
  shadeBins = 2,
  drawConnectors = FALSE,
  widthConnectors = 0.5,
  typeConnectors = 'closed',
  endsConnectors = 'first',
  lengthConnectors = unit(0.01, 'npc'),
  colConnectors = 'grey10',
  arrowheads = TRUE,
  hline = NULL,
  hlineType = 'longdash',
  hlineCol = 'black',
  hlineWidth = 0.4,
  vline = NULL,
  vlineType = 'longdash',
  vlineCol = 'black',
  vlineWidth = 0.4,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  border = 'partial',
  borderWidth = 0.8,
  borderColour = 'black', 
  raster = FALSE)
{
  if(!is.numeric(toptable[[x]])) {
    stop(paste(x, ' is not numeric!', sep=''))
  }

  if(!is.numeric(toptable[[y]])) {
    stop(paste(y, ' is not numeric!', sep=''))
  }
  
  if (raster) {
    geom_point <- geom_point_rast
  }

  i <- xvals <- yvals <- Sig <- NULL

  toptable <- as.data.frame(toptable)
  toptable$Sig <- 'NS'
  toptable$Sig[(abs(toptable[[x]]) > FCcutoff)] <- 'FC'
  toptable$Sig[(toptable[[y]] < pCutoff)] <- 'P'
  toptable$Sig[(toptable[[y]] < pCutoff) &
    (abs(toptable[[x]]) > FCcutoff)] <- 'FC_P'
  toptable$Sig <- factor(toptable$Sig,
    levels=c('NS','FC','P','FC_P'))

  # some software programs return 0 for very low p-values
  # These throw an error in EnhancedVolcano
  # Detect these, issue warning, and convert these to
  # machine-lowest value possible
  #####
  # New functionality in > v1.2:
  # Now convert to 10^-1 lower than lowest non-zero p-value
  if (min(toptable[[y]], na.rm=TRUE) == 0) {
    # <= v1.2
    #warning(paste("One or more P values is 0.",
    #  "Converting to minimum possible value..."),
    #  call. = FALSE)
    #toptable[which(toptable[[y]] == 0), y] <- .Machine$double.xmin
    warning(paste('One or more p-values is 0.',
      'Converting to 10^-1 * current',
      'lowest non-zero p-value...'),
      call. = FALSE)
    toptable[which(toptable[[y]] == 0), y] <- min(
      toptable[which(toptable[[y]] != 0), y],
      na.rm = TRUE) * 10^-1
  }

  toptable$lab <- lab
  toptable$xvals <- toptable[[x]]
  toptable$yvals <- toptable[[y]]

  # If user has supplied values in selectLab, convert labels to
  # NA and then re-set with those in selectLab
  if (!is.null(selectLab)) {
    names.new <- rep(NA, length(toptable$lab))
    indices <- which(toptable$lab %in% selectLab)
    names.new[indices] <- toptable$lab[indices]
    toptable$lab <- names.new
  }

  # create a base theme that will later be modified
  th <- theme_bw(base_size = 24) +

    theme(
      legend.background = element_rect(),

      # title, subtitle, and caption
      plot.title = element_text(
        angle = 0,
        size = titleLabSize,
        face = 'bold',
        vjust = 1),
      plot.subtitle = element_text(
        angle = 0,
        size = subtitleLabSize,
        face = 'plain',
        vjust = 1),
      plot.caption = element_text(
        angle = 0,
        size = captionLabSize,
        face = 'plain',
        vjust = 1),

      # axis text
      axis.text.x = element_text(
        angle = 0,
        size = axisLabSize,
        vjust = 1),
      axis.text.y = element_text(
        angle = 0,
        size = axisLabSize,
        vjust = 0.5),
      axis.title = element_text(
        size = axisLabSize),

      # legend
      legend.position = legendPosition,
      legend.key = element_blank(),
      legend.key.size = unit(0.5, 'cm'),
      legend.text = element_text(
        size = legendLabSize),
      title = element_text(
        size = legendLabSize),
      legend.title = element_blank())

  # Create the plot object differently based on whether colCustom 
  # and shapeCustom are NULL or not. This helps to avoid messing up
  # the legend.
  #
  # 1, both colCustom and shapeCustom are activated
  if (!is.null(colCustom) & !is.null(shapeCustom)) {

    plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

      # over-ride legend icon sizes for colour and shape.
      # guide_legends are separate for colour and shape;
      # so, legends will be drawn separate
      guides(
        colour = guide_legend(
          order = 1,
          override.aes = list(
            size = legendIconSize)),
        shape = guide_legend(
          order = 2,
          override.aes = list(
            size = legendIconSize))) +

      # include new shape and colour encodings as aes
      geom_point(
        aes(
          color = factor(names(colCustom)),
          shape = factor(names(shapeCustom))),
        alpha = colAlpha,
        size = pointSize,
        na.rm = TRUE) +

      # specify the colour and shape with the supplied encoding
      scale_color_manual(values = colCustom) +
      scale_shape_manual(values = shapeCustom)

  # 2, only colCustom is activated and 'shape' has just a single value
  } else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 1) {

    plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

      # over-ride legend icon sizes for colour and shape.
      # guide_legends are separate for colour and shape;
      # so, legends will be drawn separate IF shape is also
      # included as aes to geom_point (it is not, here)
      guides(
        colour = guide_legend(
          order = 1,
          override.aes = list(
            size = legendIconSize)),
        shape = guide_legend(
          order = 2,
          override.aes = list(
            size = legendIconSize))) +

      # include new colour encodings as aes.
      # 'shape' is included, but outside aes
      geom_point(
        aes(
          color = factor(names(colCustom))),
        alpha = colAlpha,
        shape = shape,
        size = pointSize,
        na.rm = TRUE) +

      # specify the colour with the supplied encoding
      scale_color_manual(values = colCustom) +

      # 'shape' is not included as aes. Specifying guide = TRUE
      # here will result in legends merging
      scale_shape_manual(guide = TRUE)

  # 3, only colCustom is activated and 'shape' has 4 values
  } else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 4) {

    plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

      # over-ride legend icon sizes for colour and shape.
      # guide_legends are separate for colour and shape;
      # so, legends will be drawn separate
      guides(
        colour = guide_legend(
          order = 1,
          override.aes = list(
            size = legendIconSize)),
        shape = guide_legend(
          order = 2,
          override.aes = list(
            size = legendIconSize))) +

      # include new colour encodings as aes.
      # 'shape' is included in aes and mapped to 4
      # categories of NS, FC, P, FC_P
      geom_point(
        aes(
          color = factor(names(colCustom)),
          shape = Sig),
        alpha = colAlpha,
        size = pointSize,
        na.rm = TRUE) +

      # specify the colour with the supplied encoding
      scale_color_manual(values = colCustom) +

      # as it is included as aes, a separate legend
      # for 'shape' will be drawn. Here, over-ride that
      # legend
      scale_shape_manual(
        values = c(
          NS = shape[1],
          FC = shape[2],
          P = shape[3],
          FC_P = shape[4]),
        labels = c(
          NS = legendLabels[1],
          FC = legendLabels[2],
          P = legendLabels[3],
          FC_P = legendLabels[4]),
        guide = TRUE,
        drop = legendDropLevels)

  # 4, only shapeCustom is activated
  } else if (is.null(colCustom) & !is.null(shapeCustom)) {

    if (is.null(colGradient)) {

      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        # over-ride legend icon sizes for colour and shape.
        # guide_legends are separate for colour and shape;
        # so, legends will be drawn separate
        guides(
          colour = guide_legend(
            order = 1,
            override.aes = list(
              size = legendIconSize)),
          shape = guide_legend(
            order = 2,
            override.aes = list(
              size = legendIconSize))) +

        # include new shape encodings as aes.
        # Standard colour for NS, FC, P, FC_P,
        # are added to aes, too.
        geom_point(
          aes(
            color = Sig,
            shape = factor(names(shapeCustom))),
          alpha = colAlpha,
          size = pointSize,
          na.rm = TRUE) +

        # as it is included as aes, a separate legend
        # for 'colour' will be drawn. Here, over-ride that
        # legend
        scale_color_manual(
          values = c(
            NS = col[1],
            FC = col[2],
            P = col[3],
            FC_P = col[4]),
          labels = c(
            NS = legendLabels[1],
            FC = legendLabels[2],
            P = legendLabels[3],
            FC_P = legendLabels[4]),
          drop = legendDropLevels) +

        # specify the shape with the supplied encoding
        scale_shape_manual(values = shapeCustom)

    } else {

      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        # over-ride legend icon sizes for colour and shape.
        # guide_legends are separate for colour and shape;
        # so, legends will be drawn separate
        guides(
          shape = guide_legend(
            order = 2,
            override.aes = list(
              size = legendIconSize))) +

        # include new shape encodings as aes.
        # Standard colour for NS, FC, P, FC_P,
        # are added to aes, too.
        geom_point(
          aes(
            color = Sig,
            shape = factor(names(shapeCustom))),
          alpha = colAlpha,
          size = pointSize,
          na.rm = TRUE) +

        scale_colour_gradient(
          low = colGradient[1],
          high = colGradient[2],
          limits = colGradientLimits,
          breaks = colGradientBreaks,
          labels = colGradientLabels)

        # specify the shape with the supplied encoding
        scale_shape_manual(values = shapeCustom)

    }

  # 5, both colCustom and shapeCustom are null;
  # only a single shape value specified
  } else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 1) {

    if (is.null(colGradient)) {

      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        # over-ride legend icon sizes for colour and shape.
        # including 'shape' in the colour guide_legend here
        # results in the legends merging
        guides(colour = guide_legend(
          order = 1,
          override.aes = list(
            shape = shape,
            size = legendIconSize))) +

        geom_point(
          aes(color = Sig),
          alpha = colAlpha,
          shape = shape,
          size = pointSize,
          na.rm = TRUE) +

        scale_color_manual(
          values = c(
            NS = col[1],
            FC = col[2],
            P = col[3],
            FC_P = col[4]),
          labels = c(
            NS = legendLabels[1],
            FC = legendLabels[2],
            P = legendLabels[3],
            FC_P = legendLabels[4]),
          drop = legendDropLevels)

    } else {

      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        geom_point(
          aes(color = yvals),
          alpha = colAlpha,
          shape = shape,
          size = pointSize,
          na.rm = TRUE) +

        scale_colour_gradient(
          low = colGradient[1],
          high = colGradient[2],
          limits = colGradientLimits,
          breaks = colGradientBreaks,
          labels = colGradientLabels)
    }

  # 6, both colCustom and shapeCustom are null;
  # four shape values are specified
  } else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 4) {

    if (is.null(colGradient)) {
      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        # over-ride legend icon sizes for colour and shape.
        # including 'shape' in the colour guide_legend here
        # results in the legends merging
        guides(colour = guide_legend(
          order = 1,
          override.aes = list(
            shape = c(
              NS = shape[1],
              FC = shape[2],
              P = shape[3],
              FC_P = shape[4]),
            size = legendIconSize))) +

        geom_point(
          aes(
            color = Sig,
            shape = Sig),
          alpha = colAlpha,
          size = pointSize,
          na.rm = TRUE) +

        scale_color_manual(
          values = c(
            NS = col[1],
            FC = col[2],
            P = col[3],
            FC_P = col[4]),
          labels = c(
            NS = legendLabels[1],
            FC = legendLabels[2],
            P = legendLabels[3],
            FC_P = legendLabels[4]),
          drop = legendDropLevels) +

        scale_shape_manual(
          values = c(
            NS = shape[1],
            FC = shape[2],
            P = shape[3],
            FC_P = shape[4]),
          guide = FALSE,
          drop = legendDropLevels)

    } else {

      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        geom_point(
          aes(
            color = yvals,
            shape = Sig),
          alpha = colAlpha,
          size = pointSize,
          na.rm = TRUE) +

        scale_colour_gradient(
          low = colGradient[1],
          high = colGradient[2],
          limits = colGradientLimits,
          breaks = colGradientBreaks,
          labels = colGradientLabels) +

        scale_shape_manual(
          values = c(
            NS = shape[1],
            FC = shape[2],
            P = shape[3],
            FC_P = shape[4]),
          guide = FALSE,
          drop = legendDropLevels)

    }
  }

  # add more elements to the plot
  plot <- plot +

    xlab(xlab) +
    ylab(ylab) +

    xlim(xlim[1], xlim[2]) +
    ylim(ylim[1], ylim[2]) +

    geom_vline(xintercept = c(-FCcutoff, FCcutoff),
      linetype = cutoffLineType,
      colour = cutoffLineCol,
      size = cutoffLineWidth) +

    geom_hline(yintercept = -log10(pCutoff),
      linetype = cutoffLineType,
      colour = cutoffLineCol,
      size = cutoffLineWidth)

  # add elements to the plot for title, subtitle, caption
  plot <- plot + labs(title = title, 
    subtitle = subtitle, caption = caption)

  # add elements to the plot for vlines and hlines
  if (!is.null(vline)) {
    plot <- plot + geom_vline(xintercept = vline,
      linetype = vlineType,
      colour = vlineCol,
      size = vlineWidth)
  }
  if (!is.null(hline)) {
    plot <- plot + geom_hline(yintercept = -log10(hline),
      linetype = hlineType,
      colour = hlineCol,
      size = hlineWidth)
  }

  # Border around plot
  if (border == 'full') {
    plot <- plot + theme(panel.border = element_rect(
      colour = borderColour, fill = NA, size = borderWidth))
  } else if (border == 'partial') {
    plot <- plot + theme(axis.line = element_line(
      size = borderWidth, colour = borderColour),
      panel.border = element_blank(),
      panel.background = element_blank())
  } else {
    stop('Unrecognised value passed to \'border\'. Must be \'full\' or \'partial\'')
  }

  # Gridlines
  if (gridlines.major) {
    plot <- plot + theme(panel.grid.major = element_line())
  } else {
    plot <- plot + theme(panel.grid.major = element_blank())
  }
  if (gridlines.minor) {
    plot <- plot + theme(panel.grid.minor = element_line())
  } else {
    plot <- plot + theme(panel.grid.minor = element_blank())
  }

  # user has specified to draw with geom_text or geom_label?
  if (!boxedLabels) {

    # For labeling with geom_[text|label]_repel and
    # geom_[text|label] with check_overlap = TRUE, 4 possible
    # scenarios can arise
    if (drawConnectors && is.null(selectLab)) {

      if (arrowheads) {
        arr <- arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors)
      } else {
        arr <- NULL
      }

      plot <- plot + geom_text_repel(
        data = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff),
        aes(label = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff)[["lab"]]),
        size = labSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arr,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    } else if (drawConnectors && !is.null(selectLab)) {

      if (arrowheads) {
        arr <- arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors)
      } else {
        arr <- NULL
      }

      plot <- plot + geom_text_repel(
        data = subset(toptable,
          !is.na(toptable[['lab']])),
        aes(label = subset(toptable,
          !is.na(toptable[['lab']]))[['lab']]),
        size = labSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arr,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    } else if (!drawConnectors && !is.null(selectLab)) {

      plot <- plot + geom_text(
        data = subset(toptable,
          !is.na(toptable[['lab']])),
        aes(
          label = subset(toptable,
            !is.na(toptable[['lab']]))[['lab']]),
        size = labSize,
        check_overlap = TRUE,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    } else if (!drawConnectors && is.null(selectLab)) {

      plot <- plot + geom_text(
        data = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff),
        aes(label = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff)[['lab']]),
        size = labSize,
        check_overlap = TRUE,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)
    }

  } else {

    # For labeling with geom_[text|label]_repel and
    # geom_[text|label] with check_overlap = TRUE, 4 possible
    # scenarios can arise
    if (drawConnectors && is.null(selectLab)) {

      if (arrowheads) {
        arr <- arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors)
      } else {
        arr <- NULL
      }

      plot <- plot + geom_label_repel(
        data = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff),
        aes(label = subset(toptable,
          toptable[[y]]<pCutoff &
            abs(toptable[[x]]) > FCcutoff)[['lab']]),
        size = labSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arr,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    } else if (drawConnectors && !is.null(selectLab)) {

      if (arrowheads) {
        arr <- arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors)
      } else {
        arr <- NULL
      }

      plot <- plot + geom_label_repel(
        data = subset(toptable,
          !is.na(toptable[['lab']])),
        aes(label = subset(toptable,
          !is.na(toptable[['lab']]))[['lab']]),
        size = labSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arr,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    } else if (!drawConnectors && !is.null(selectLab)) {

      plot <- plot + geom_label(
        data = subset(toptable,
          !is.na(toptable[["lab"]])),
        aes(
          label = subset(toptable,
            !is.na(toptable[['lab']]))[['lab']]),
        size = labSize,
        #check_overlap = TRUE,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    } else if (!drawConnectors && is.null(selectLab)) {

      plot <- plot + geom_label(
        data = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff),
        aes(label = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff)[['lab']]),
        size = labSize,
        #check_overlap = TRUE,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    }
  }

  # encircle
  if (!is.null(encircle)) {
    plot <- plot + 
      geom_encircle(
        data = subset(toptable,
          rownames(toptable) %in% encircle),
        colour = encircleCol,
        fill = encircleFill,
        alpha = encircleAlpha,
        size = encircleSize,
        show.legend = FALSE,
        na.rm = TRUE)
  }

  # shade
  if (!is.null(shade)) {
    plot <- plot + 
      stat_density2d(
        data = subset(toptable,
          rownames(toptable) %in% shade),
        fill = shadeFill,
        alpha = shadeAlpha,
        geom = 'polygon',
        contour = TRUE,
        size = shadeSize,
        bins = shadeBins,
        show.legend = FALSE,
        na.rm = TRUE)
  }

  return(plot)
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.a6be4b2e-17f5-488e-8c25-efe361a38e0a"),
    naive_deg_analysis=Input(rid="ri.foundry.main.dataset.862f8327-dcb8-414e-b8ae-fb74aaa3a15b")
)
volcano_d2vd1 <- function(naive_deg_analysis) {
    # image: png

    # Changelog
    # 2020-10-29 Add support for pval == 0

    suppressMessages(library(stringr))
    suppressMessages(library(tidyverse))
    suppressMessages(library(ggplot2))
    
    #Class handling
    df <- naive_deg_analysis
    if(class(df) == "SparkDataFrame"){
        df <- SparkR::collect(df)
    }
  
    # Log fold-change column
    # p-value column
    # Position of p-value threshold line #default 0.05
    # Position of log fold-change threshold lines #default 1 lfc
    # Label column #default gene names
    # Specific genes to label
    # Number of genes to label
    
    #label topN genes or specific genes
    
    label.col <- "Gene"
    sig.col <- "d2-d1_adjpval"
    lfc.col <- "d2-d1_logFC"
    columns_of_interest <- c(label.col,lfc.col,sig.col)

    df <- df %>% dplyr::select(one_of(columns_of_interest)) %>% dplyr::filter(!is.na(!!rlang::sym(lfc.col))) #seurat introduces NAs

    change_lfc_name <- "log2FC"
    change_sig_name <- "FDR"
    colnames(df) <- c(label.col,change_lfc_name,change_sig_name)

    cat(paste0("Genes in initial dataset: ", nrow(df),"\n"))

    #Select top genes by logFC or Significance
    no_genes_to_label <- 30
    value_to_sort_the_output_dataset <- "p-value"
    if (value_to_sort_the_output_dataset=="fold-change") {
        df %>% dplyr::arrange(desc(abs(!!rlang::sym(change_lfc_name)))) -> df
    } else if (value_to_sort_the_output_dataset=="p-value") {
        df %>% dplyr::arrange(!!rlang::sym(change_sig_name)) -> df
    }
    genes_to_label <- as.character(df[1:no_genes_to_label,label.col])

    additional_labels <- ""
    additional_labels <- unlist(str_split(additional_labels,","))
    filter <- additional_labels %in% df[,label.col]
    additional_labels <- additional_labels[filter]
    missing_labels <- additional_labels[!filter]

    if(!is.na(missing_labels)){
        cat("Could not find:\n")
        print(missing_labels)
    }
    use_only_addition_labels <- FALSE
    if(use_only_addition_labels){
        genes_to_label <- additional_labels
    }else{
        genes_to_label <- unique(append(genes_to_label, additional_labels))
    }

    # Plot variables
    title <- "Volcano Plots"
    subtitle <- "Vaccination 1 ( 1 dpv vs pre)"

    pCutoff  = 0.05
    FCcutoff = 1.0

    significant=as.vector(table(abs( df[,change_lfc_name] ) > FCcutoff &
                                     df[,change_sig_name]   < pCutoff))[2]

    # fix pvalue == 0
    shapeCustom <- rep(19,nrow(df))
    maxy <-  max(-log10(df[[change_sig_name]]), na.rm=TRUE)
    if(0 > 0){
        maxy <- 0
    }
    
    cat(paste0("Maxy: ",maxy,"\n"))
    if(maxy == Inf){
        # Sometimes, pvalues == 0
        keep <- df[[change_sig_name]] > 0
        df[[change_sig_name]][!keep] <- min(df[[change_sig_name]][keep])
        shapeCustom[!keep] <- 17

        maxy <- -log10(min(df[[change_sig_name]][keep]))
        cat("Some p-values equal zero. Adjusting y-limits.\n")
        cat(paste0("Maxy adjusted: ",maxy,"\n"))

    }

    # By default, nothing will be greater than maxy. User can set this value lower
    keep <- -log10(df[[change_sig_name]]) <= maxy
    df[[change_sig_name]][!keep] <- maxy
    shapeCustom[!keep] <- 17

    names(shapeCustom)<- rep("Exact",length(shapeCustom))
    names(shapeCustom)[shapeCustom == 17] <- "Adjusted"
    
    #Remove if nothin' doin'
    if(all(shapeCustom == 19)){
        shapeCustom <- NULL
    }
    
    maxy <- ceiling(maxy)

    use_custom_xlab <- FALSE
    if(use_custom_xlab){
        xlab <- gsub("_"," ",change_lfc_name)
    }else{
        xlab <- bquote(~Log[2]~ "fold change")
    }

    xlim_additional <- 0
    ylim_additional <- 0
    axisLabSize <- 24
    labSize <- 4
    pointSize <- 2

    imageWidth = 3000
    imageHeight = 3000
    dpi = 300

    png(
      filename=graphicsFile,
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")

    p <- EnhancedVolcano(   df,x=change_lfc_name,y=change_sig_name,
                            lab=df[,label.col],
                            selectLab = genes_to_label,
                            title=paste0(title," (Significant=",significant,")"),
                            subtitle <- subtitle,
                            xlab=xlab,
                            ylab=bquote(~-Log[10]~.(change_sig_name)),
                            xlim=c(floor(min(df[,change_lfc_name])) - xlim_additional,ceiling(max(df[,change_lfc_name]))+ xlim_additional),
                            ylim=c(0, maxy + ylim_additional),
                            pCutoff=pCutoff,
                            FCcutoff=FCcutoff,
                            axisLabSize=axisLabSize,
                            labSize=labSize,
                            pointSize=pointSize,
                            shapeCustom=shapeCustom
                            )
    print(p)

    df$rank <- -log10(df[[change_sig_name]]) * sign(df[[change_lfc_name]]) #sig already -log10
    df <- df %>% arrange(desc(rank))
    return(df)
}

# This is directly copied from https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html v1.6.0
EnhancedVolcano <- function(
  toptable,
  lab,
  x,
  y,
  selectLab = NULL,
  xlim = c(min(toptable[[x]], na.rm=TRUE) - 1.5,
    max(toptable[[x]], na.rm=TRUE) + 1.5),
  ylim = c(0, max(-log10(toptable[[y]]), na.rm=TRUE) + 5),
  xlab = bquote(~Log[2]~ "fold change"),
  ylab = bquote(~-Log[10]~italic(P)),
  axisLabSize = 18,
  title = 'Volcano plot',
  subtitle = bquote(italic(EnhancedVolcano)),
  caption = paste0('total = ', nrow(toptable), ' variables'),
  titleLabSize = 18,
  subtitleLabSize = 14,
  captionLabSize = 14,
  pCutoff = 10e-6,
  FCcutoff = 1.0,
  cutoffLineType = 'longdash',
  cutoffLineCol = 'black',
  cutoffLineWidth = 0.4,
  pointSize = 2.0,
  labSize = 5.0,
  labCol = 'black',
  labFace = 'plain',
  labhjust = 0.5,
  labvjust = 1.5,
  boxedLabels = FALSE,
  shape = 19,
  shapeCustom = NULL,
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  colCustom = NULL,
  colAlpha = 1/2,
  colGradient = NULL,
  colGradientBreaks = c(pCutoff, 1.0),
  colGradientLabels = c('0', '1.0'),
  colGradientLimits = c(0, 1.0),
  legendLabels = c('NS', expression(Log[2]~FC),
    'p-value', expression(p-value~and~log[2]~FC)),
  legendPosition = 'top',
  legendLabSize = 14,
  legendIconSize = 5.0,
  legendDropLevels = TRUE,
  encircle = NULL,
  encircleCol = 'black',
  encircleFill = 'pink',
  encircleAlpha = 3/4,
  encircleSize = 2.5,
  shade = NULL,
  shadeFill = 'grey',
  shadeAlpha = 1/2,
  shadeSize = 0.01,
  shadeBins = 2,
  drawConnectors = FALSE,
  widthConnectors = 0.5,
  typeConnectors = 'closed',
  endsConnectors = 'first',
  lengthConnectors = unit(0.01, 'npc'),
  colConnectors = 'grey10',
  arrowheads = TRUE,
  hline = NULL,
  hlineType = 'longdash',
  hlineCol = 'black',
  hlineWidth = 0.4,
  vline = NULL,
  vlineType = 'longdash',
  vlineCol = 'black',
  vlineWidth = 0.4,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  border = 'partial',
  borderWidth = 0.8,
  borderColour = 'black', 
  raster = FALSE)
{
  if(!is.numeric(toptable[[x]])) {
    stop(paste(x, ' is not numeric!', sep=''))
  }

  if(!is.numeric(toptable[[y]])) {
    stop(paste(y, ' is not numeric!', sep=''))
  }
  
  if (raster) {
    geom_point <- geom_point_rast
  }

  i <- xvals <- yvals <- Sig <- NULL

  toptable <- as.data.frame(toptable)
  toptable$Sig <- 'NS'
  toptable$Sig[(abs(toptable[[x]]) > FCcutoff)] <- 'FC'
  toptable$Sig[(toptable[[y]] < pCutoff)] <- 'P'
  toptable$Sig[(toptable[[y]] < pCutoff) &
    (abs(toptable[[x]]) > FCcutoff)] <- 'FC_P'
  toptable$Sig <- factor(toptable$Sig,
    levels=c('NS','FC','P','FC_P'))

  # some software programs return 0 for very low p-values
  # These throw an error in EnhancedVolcano
  # Detect these, issue warning, and convert these to
  # machine-lowest value possible
  #####
  # New functionality in > v1.2:
  # Now convert to 10^-1 lower than lowest non-zero p-value
  if (min(toptable[[y]], na.rm=TRUE) == 0) {
    # <= v1.2
    #warning(paste("One or more P values is 0.",
    #  "Converting to minimum possible value..."),
    #  call. = FALSE)
    #toptable[which(toptable[[y]] == 0), y] <- .Machine$double.xmin
    warning(paste('One or more p-values is 0.',
      'Converting to 10^-1 * current',
      'lowest non-zero p-value...'),
      call. = FALSE)
    toptable[which(toptable[[y]] == 0), y] <- min(
      toptable[which(toptable[[y]] != 0), y],
      na.rm = TRUE) * 10^-1
  }

  toptable$lab <- lab
  toptable$xvals <- toptable[[x]]
  toptable$yvals <- toptable[[y]]

  # If user has supplied values in selectLab, convert labels to
  # NA and then re-set with those in selectLab
  if (!is.null(selectLab)) {
    names.new <- rep(NA, length(toptable$lab))
    indices <- which(toptable$lab %in% selectLab)
    names.new[indices] <- toptable$lab[indices]
    toptable$lab <- names.new
  }

  # create a base theme that will later be modified
  th <- theme_bw(base_size = 24) +

    theme(
      legend.background = element_rect(),

      # title, subtitle, and caption
      plot.title = element_text(
        angle = 0,
        size = titleLabSize,
        face = 'bold',
        vjust = 1),
      plot.subtitle = element_text(
        angle = 0,
        size = subtitleLabSize,
        face = 'plain',
        vjust = 1),
      plot.caption = element_text(
        angle = 0,
        size = captionLabSize,
        face = 'plain',
        vjust = 1),

      # axis text
      axis.text.x = element_text(
        angle = 0,
        size = axisLabSize,
        vjust = 1),
      axis.text.y = element_text(
        angle = 0,
        size = axisLabSize,
        vjust = 0.5),
      axis.title = element_text(
        size = axisLabSize),

      # legend
      legend.position = legendPosition,
      legend.key = element_blank(),
      legend.key.size = unit(0.5, 'cm'),
      legend.text = element_text(
        size = legendLabSize),
      title = element_text(
        size = legendLabSize),
      legend.title = element_blank())

  # Create the plot object differently based on whether colCustom 
  # and shapeCustom are NULL or not. This helps to avoid messing up
  # the legend.
  #
  # 1, both colCustom and shapeCustom are activated
  if (!is.null(colCustom) & !is.null(shapeCustom)) {

    plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

      # over-ride legend icon sizes for colour and shape.
      # guide_legends are separate for colour and shape;
      # so, legends will be drawn separate
      guides(
        colour = guide_legend(
          order = 1,
          override.aes = list(
            size = legendIconSize)),
        shape = guide_legend(
          order = 2,
          override.aes = list(
            size = legendIconSize))) +

      # include new shape and colour encodings as aes
      geom_point(
        aes(
          color = factor(names(colCustom)),
          shape = factor(names(shapeCustom))),
        alpha = colAlpha,
        size = pointSize,
        na.rm = TRUE) +

      # specify the colour and shape with the supplied encoding
      scale_color_manual(values = colCustom) +
      scale_shape_manual(values = shapeCustom)

  # 2, only colCustom is activated and 'shape' has just a single value
  } else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 1) {

    plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

      # over-ride legend icon sizes for colour and shape.
      # guide_legends are separate for colour and shape;
      # so, legends will be drawn separate IF shape is also
      # included as aes to geom_point (it is not, here)
      guides(
        colour = guide_legend(
          order = 1,
          override.aes = list(
            size = legendIconSize)),
        shape = guide_legend(
          order = 2,
          override.aes = list(
            size = legendIconSize))) +

      # include new colour encodings as aes.
      # 'shape' is included, but outside aes
      geom_point(
        aes(
          color = factor(names(colCustom))),
        alpha = colAlpha,
        shape = shape,
        size = pointSize,
        na.rm = TRUE) +

      # specify the colour with the supplied encoding
      scale_color_manual(values = colCustom) +

      # 'shape' is not included as aes. Specifying guide = TRUE
      # here will result in legends merging
      scale_shape_manual(guide = TRUE)

  # 3, only colCustom is activated and 'shape' has 4 values
  } else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 4) {

    plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

      # over-ride legend icon sizes for colour and shape.
      # guide_legends are separate for colour and shape;
      # so, legends will be drawn separate
      guides(
        colour = guide_legend(
          order = 1,
          override.aes = list(
            size = legendIconSize)),
        shape = guide_legend(
          order = 2,
          override.aes = list(
            size = legendIconSize))) +

      # include new colour encodings as aes.
      # 'shape' is included in aes and mapped to 4
      # categories of NS, FC, P, FC_P
      geom_point(
        aes(
          color = factor(names(colCustom)),
          shape = Sig),
        alpha = colAlpha,
        size = pointSize,
        na.rm = TRUE) +

      # specify the colour with the supplied encoding
      scale_color_manual(values = colCustom) +

      # as it is included as aes, a separate legend
      # for 'shape' will be drawn. Here, over-ride that
      # legend
      scale_shape_manual(
        values = c(
          NS = shape[1],
          FC = shape[2],
          P = shape[3],
          FC_P = shape[4]),
        labels = c(
          NS = legendLabels[1],
          FC = legendLabels[2],
          P = legendLabels[3],
          FC_P = legendLabels[4]),
        guide = TRUE,
        drop = legendDropLevels)

  # 4, only shapeCustom is activated
  } else if (is.null(colCustom) & !is.null(shapeCustom)) {

    if (is.null(colGradient)) {

      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        # over-ride legend icon sizes for colour and shape.
        # guide_legends are separate for colour and shape;
        # so, legends will be drawn separate
        guides(
          colour = guide_legend(
            order = 1,
            override.aes = list(
              size = legendIconSize)),
          shape = guide_legend(
            order = 2,
            override.aes = list(
              size = legendIconSize))) +

        # include new shape encodings as aes.
        # Standard colour for NS, FC, P, FC_P,
        # are added to aes, too.
        geom_point(
          aes(
            color = Sig,
            shape = factor(names(shapeCustom))),
          alpha = colAlpha,
          size = pointSize,
          na.rm = TRUE) +

        # as it is included as aes, a separate legend
        # for 'colour' will be drawn. Here, over-ride that
        # legend
        scale_color_manual(
          values = c(
            NS = col[1],
            FC = col[2],
            P = col[3],
            FC_P = col[4]),
          labels = c(
            NS = legendLabels[1],
            FC = legendLabels[2],
            P = legendLabels[3],
            FC_P = legendLabels[4]),
          drop = legendDropLevels) +

        # specify the shape with the supplied encoding
        scale_shape_manual(values = shapeCustom)

    } else {

      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        # over-ride legend icon sizes for colour and shape.
        # guide_legends are separate for colour and shape;
        # so, legends will be drawn separate
        guides(
          shape = guide_legend(
            order = 2,
            override.aes = list(
              size = legendIconSize))) +

        # include new shape encodings as aes.
        # Standard colour for NS, FC, P, FC_P,
        # are added to aes, too.
        geom_point(
          aes(
            color = Sig,
            shape = factor(names(shapeCustom))),
          alpha = colAlpha,
          size = pointSize,
          na.rm = TRUE) +

        scale_colour_gradient(
          low = colGradient[1],
          high = colGradient[2],
          limits = colGradientLimits,
          breaks = colGradientBreaks,
          labels = colGradientLabels)

        # specify the shape with the supplied encoding
        scale_shape_manual(values = shapeCustom)

    }

  # 5, both colCustom and shapeCustom are null;
  # only a single shape value specified
  } else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 1) {

    if (is.null(colGradient)) {

      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        # over-ride legend icon sizes for colour and shape.
        # including 'shape' in the colour guide_legend here
        # results in the legends merging
        guides(colour = guide_legend(
          order = 1,
          override.aes = list(
            shape = shape,
            size = legendIconSize))) +

        geom_point(
          aes(color = Sig),
          alpha = colAlpha,
          shape = shape,
          size = pointSize,
          na.rm = TRUE) +

        scale_color_manual(
          values = c(
            NS = col[1],
            FC = col[2],
            P = col[3],
            FC_P = col[4]),
          labels = c(
            NS = legendLabels[1],
            FC = legendLabels[2],
            P = legendLabels[3],
            FC_P = legendLabels[4]),
          drop = legendDropLevels)

    } else {

      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        geom_point(
          aes(color = yvals),
          alpha = colAlpha,
          shape = shape,
          size = pointSize,
          na.rm = TRUE) +

        scale_colour_gradient(
          low = colGradient[1],
          high = colGradient[2],
          limits = colGradientLimits,
          breaks = colGradientBreaks,
          labels = colGradientLabels)
    }

  # 6, both colCustom and shapeCustom are null;
  # four shape values are specified
  } else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 4) {

    if (is.null(colGradient)) {
      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        # over-ride legend icon sizes for colour and shape.
        # including 'shape' in the colour guide_legend here
        # results in the legends merging
        guides(colour = guide_legend(
          order = 1,
          override.aes = list(
            shape = c(
              NS = shape[1],
              FC = shape[2],
              P = shape[3],
              FC_P = shape[4]),
            size = legendIconSize))) +

        geom_point(
          aes(
            color = Sig,
            shape = Sig),
          alpha = colAlpha,
          size = pointSize,
          na.rm = TRUE) +

        scale_color_manual(
          values = c(
            NS = col[1],
            FC = col[2],
            P = col[3],
            FC_P = col[4]),
          labels = c(
            NS = legendLabels[1],
            FC = legendLabels[2],
            P = legendLabels[3],
            FC_P = legendLabels[4]),
          drop = legendDropLevels) +

        scale_shape_manual(
          values = c(
            NS = shape[1],
            FC = shape[2],
            P = shape[3],
            FC_P = shape[4]),
          guide = FALSE,
          drop = legendDropLevels)

    } else {

      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        geom_point(
          aes(
            color = yvals,
            shape = Sig),
          alpha = colAlpha,
          size = pointSize,
          na.rm = TRUE) +

        scale_colour_gradient(
          low = colGradient[1],
          high = colGradient[2],
          limits = colGradientLimits,
          breaks = colGradientBreaks,
          labels = colGradientLabels) +

        scale_shape_manual(
          values = c(
            NS = shape[1],
            FC = shape[2],
            P = shape[3],
            FC_P = shape[4]),
          guide = FALSE,
          drop = legendDropLevels)

    }
  }

  # add more elements to the plot
  plot <- plot +

    xlab(xlab) +
    ylab(ylab) +

    xlim(xlim[1], xlim[2]) +
    ylim(ylim[1], ylim[2]) +

    geom_vline(xintercept = c(-FCcutoff, FCcutoff),
      linetype = cutoffLineType,
      colour = cutoffLineCol,
      size = cutoffLineWidth) +

    geom_hline(yintercept = -log10(pCutoff),
      linetype = cutoffLineType,
      colour = cutoffLineCol,
      size = cutoffLineWidth)

  # add elements to the plot for title, subtitle, caption
  plot <- plot + labs(title = title, 
    subtitle = subtitle, caption = caption)

  # add elements to the plot for vlines and hlines
  if (!is.null(vline)) {
    plot <- plot + geom_vline(xintercept = vline,
      linetype = vlineType,
      colour = vlineCol,
      size = vlineWidth)
  }
  if (!is.null(hline)) {
    plot <- plot + geom_hline(yintercept = -log10(hline),
      linetype = hlineType,
      colour = hlineCol,
      size = hlineWidth)
  }

  # Border around plot
  if (border == 'full') {
    plot <- plot + theme(panel.border = element_rect(
      colour = borderColour, fill = NA, size = borderWidth))
  } else if (border == 'partial') {
    plot <- plot + theme(axis.line = element_line(
      size = borderWidth, colour = borderColour),
      panel.border = element_blank(),
      panel.background = element_blank())
  } else {
    stop('Unrecognised value passed to \'border\'. Must be \'full\' or \'partial\'')
  }

  # Gridlines
  if (gridlines.major) {
    plot <- plot + theme(panel.grid.major = element_line())
  } else {
    plot <- plot + theme(panel.grid.major = element_blank())
  }
  if (gridlines.minor) {
    plot <- plot + theme(panel.grid.minor = element_line())
  } else {
    plot <- plot + theme(panel.grid.minor = element_blank())
  }

  # user has specified to draw with geom_text or geom_label?
  if (!boxedLabels) {

    # For labeling with geom_[text|label]_repel and
    # geom_[text|label] with check_overlap = TRUE, 4 possible
    # scenarios can arise
    if (drawConnectors && is.null(selectLab)) {

      if (arrowheads) {
        arr <- arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors)
      } else {
        arr <- NULL
      }

      plot <- plot + geom_text_repel(
        data = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff),
        aes(label = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff)[["lab"]]),
        size = labSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arr,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    } else if (drawConnectors && !is.null(selectLab)) {

      if (arrowheads) {
        arr <- arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors)
      } else {
        arr <- NULL
      }

      plot <- plot + geom_text_repel(
        data = subset(toptable,
          !is.na(toptable[['lab']])),
        aes(label = subset(toptable,
          !is.na(toptable[['lab']]))[['lab']]),
        size = labSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arr,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    } else if (!drawConnectors && !is.null(selectLab)) {

      plot <- plot + geom_text(
        data = subset(toptable,
          !is.na(toptable[['lab']])),
        aes(
          label = subset(toptable,
            !is.na(toptable[['lab']]))[['lab']]),
        size = labSize,
        check_overlap = TRUE,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    } else if (!drawConnectors && is.null(selectLab)) {

      plot <- plot + geom_text(
        data = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff),
        aes(label = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff)[['lab']]),
        size = labSize,
        check_overlap = TRUE,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)
    }

  } else {

    # For labeling with geom_[text|label]_repel and
    # geom_[text|label] with check_overlap = TRUE, 4 possible
    # scenarios can arise
    if (drawConnectors && is.null(selectLab)) {

      if (arrowheads) {
        arr <- arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors)
      } else {
        arr <- NULL
      }

      plot <- plot + geom_label_repel(
        data = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff),
        aes(label = subset(toptable,
          toptable[[y]]<pCutoff &
            abs(toptable[[x]]) > FCcutoff)[['lab']]),
        size = labSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arr,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    } else if (drawConnectors && !is.null(selectLab)) {

      if (arrowheads) {
        arr <- arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors)
      } else {
        arr <- NULL
      }

      plot <- plot + geom_label_repel(
        data = subset(toptable,
          !is.na(toptable[['lab']])),
        aes(label = subset(toptable,
          !is.na(toptable[['lab']]))[['lab']]),
        size = labSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arr,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    } else if (!drawConnectors && !is.null(selectLab)) {

      plot <- plot + geom_label(
        data = subset(toptable,
          !is.na(toptable[["lab"]])),
        aes(
          label = subset(toptable,
            !is.na(toptable[['lab']]))[['lab']]),
        size = labSize,
        #check_overlap = TRUE,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    } else if (!drawConnectors && is.null(selectLab)) {

      plot <- plot + geom_label(
        data = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff),
        aes(label = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff)[['lab']]),
        size = labSize,
        #check_overlap = TRUE,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    }
  }

  # encircle
  if (!is.null(encircle)) {
    plot <- plot + 
      geom_encircle(
        data = subset(toptable,
          rownames(toptable) %in% encircle),
        colour = encircleCol,
        fill = encircleFill,
        alpha = encircleAlpha,
        size = encircleSize,
        show.legend = FALSE,
        na.rm = TRUE)
  }

  # shade
  if (!is.null(shade)) {
    plot <- plot + 
      stat_density2d(
        data = subset(toptable,
          rownames(toptable) %in% shade),
        fill = shadeFill,
        alpha = shadeAlpha,
        geom = 'polygon',
        contour = TRUE,
        size = shadeSize,
        bins = shadeBins,
        show.legend = FALSE,
        na.rm = TRUE)
  }

  return(plot)
}

#################################################
## Global imports and functions included below ##
#################################################

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.38ec61d7-cd64-4207-b55b-442dab1ccc8c"),
    naive_deg_analysis=Input(rid="ri.foundry.main.dataset.862f8327-dcb8-414e-b8ae-fb74aaa3a15b")
)
volcano_d8vd1 <- function(naive_deg_analysis) {
    # image: png

    # Changelog
    # 2020-10-29 Add support for pval == 0

    suppressMessages(library(stringr))
    suppressMessages(library(tidyverse))
    suppressMessages(library(ggplot2))
    
    #Class handling
    df <- naive_deg_analysis
    if(class(df) == "SparkDataFrame"){
        df <- SparkR::collect(df)
    }
  
    # Log fold-change column
    # p-value column
    # Position of p-value threshold line #default 0.05
    # Position of log fold-change threshold lines #default 1 lfc
    # Label column #default gene names
    # Specific genes to label
    # Number of genes to label
    
    #label topN genes or specific genes
    
    label.col <- "Gene"
    sig.col <- "d8-d1_adjpval"
    lfc.col <- "d8-d1_logFC"
    columns_of_interest <- c(label.col,lfc.col,sig.col)

    df <- df %>% dplyr::select(one_of(columns_of_interest)) %>% dplyr::filter(!is.na(!!rlang::sym(lfc.col))) #seurat introduces NAs

    change_lfc_name <- "log2FC"
    change_sig_name <- "FDR"
    colnames(df) <- c(label.col,change_lfc_name,change_sig_name)

    cat(paste0("Genes in initial dataset: ", nrow(df),"\n"))

    #Select top genes by logFC or Significance
    no_genes_to_label <- 30
    value_to_sort_the_output_dataset <- "p-value"
    if (value_to_sort_the_output_dataset=="fold-change") {
        df %>% dplyr::arrange(desc(abs(!!rlang::sym(change_lfc_name)))) -> df
    } else if (value_to_sort_the_output_dataset=="p-value") {
        df %>% dplyr::arrange(!!rlang::sym(change_sig_name)) -> df
    }
    genes_to_label <- as.character(df[1:no_genes_to_label,label.col])

    additional_labels <- ""
    additional_labels <- unlist(str_split(additional_labels,","))
    filter <- additional_labels %in% df[,label.col]
    additional_labels <- additional_labels[filter]
    missing_labels <- additional_labels[!filter]

    if(!is.na(missing_labels)){
        cat("Could not find:\n")
        print(missing_labels)
    }
    use_only_addition_labels <- FALSE
    if(use_only_addition_labels){
        genes_to_label <- additional_labels
    }else{
        genes_to_label <- unique(append(genes_to_label, additional_labels))
    }

    # Plot variables
    title <- "Volcano Plots"
    subtitle <- "Vaccination 1 ( 7 dpv vs pre)"

    pCutoff  = 0.05
    FCcutoff = 1.0

    significant=as.vector(table(abs( df[,change_lfc_name] ) > FCcutoff &
                                     df[,change_sig_name]   < pCutoff))[2]

    # fix pvalue == 0
    shapeCustom <- rep(19,nrow(df))
    maxy <-  max(-log10(df[[change_sig_name]]), na.rm=TRUE)
    if(0 > 0){
        maxy <- 0
    }
    
    cat(paste0("Maxy: ",maxy,"\n"))
    if(maxy == Inf){
        # Sometimes, pvalues == 0
        keep <- df[[change_sig_name]] > 0
        df[[change_sig_name]][!keep] <- min(df[[change_sig_name]][keep])
        shapeCustom[!keep] <- 17

        maxy <- -log10(min(df[[change_sig_name]][keep]))
        cat("Some p-values equal zero. Adjusting y-limits.\n")
        cat(paste0("Maxy adjusted: ",maxy,"\n"))

    }

    # By default, nothing will be greater than maxy. User can set this value lower
    keep <- -log10(df[[change_sig_name]]) <= maxy
    df[[change_sig_name]][!keep] <- maxy
    shapeCustom[!keep] <- 17

    names(shapeCustom)<- rep("Exact",length(shapeCustom))
    names(shapeCustom)[shapeCustom == 17] <- "Adjusted"
    
    #Remove if nothin' doin'
    if(all(shapeCustom == 19)){
        shapeCustom <- NULL
    }
    
    maxy <- ceiling(maxy)

    use_custom_xlab <- FALSE
    if(use_custom_xlab){
        xlab <- gsub("_"," ",change_lfc_name)
    }else{
        xlab <- bquote(~Log[2]~ "fold change")
    }

    xlim_additional <- 0
    ylim_additional <- 0
    axisLabSize <- 24
    labSize <- 4
    pointSize <- 2

    imageWidth = 3000
    imageHeight = 3000
    dpi = 300

    png(
      filename=graphicsFile,
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")

    p <- EnhancedVolcano(   df,x=change_lfc_name,y=change_sig_name,
                            lab=df[,label.col],
                            selectLab = genes_to_label,
                            title=paste0(title," (Significant=",significant,")"),
                            subtitle <- subtitle,
                            xlab=xlab,
                            ylab=bquote(~-Log[10]~.(change_sig_name)),
                            xlim=c(floor(min(df[,change_lfc_name])) - xlim_additional,ceiling(max(df[,change_lfc_name]))+ xlim_additional),
                            ylim=c(0, maxy + ylim_additional),
                            pCutoff=pCutoff,
                            FCcutoff=FCcutoff,
                            axisLabSize=axisLabSize,
                            labSize=labSize,
                            pointSize=pointSize,
                            shapeCustom=shapeCustom
                            )
    print(p)

    df$rank <- -log10(df[[change_sig_name]]) * sign(df[[change_lfc_name]]) #sig already -log10
    df <- df %>% arrange(desc(rank))
    return(df)
}

# This is directly copied from https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html v1.6.0
EnhancedVolcano <- function(
  toptable,
  lab,
  x,
  y,
  selectLab = NULL,
  xlim = c(min(toptable[[x]], na.rm=TRUE) - 1.5,
    max(toptable[[x]], na.rm=TRUE) + 1.5),
  ylim = c(0, max(-log10(toptable[[y]]), na.rm=TRUE) + 5),
  xlab = bquote(~Log[2]~ "fold change"),
  ylab = bquote(~-Log[10]~italic(P)),
  axisLabSize = 18,
  title = 'Volcano plot',
  subtitle = bquote(italic(EnhancedVolcano)),
  caption = paste0('total = ', nrow(toptable), ' variables'),
  titleLabSize = 18,
  subtitleLabSize = 14,
  captionLabSize = 14,
  pCutoff = 10e-6,
  FCcutoff = 1.0,
  cutoffLineType = 'longdash',
  cutoffLineCol = 'black',
  cutoffLineWidth = 0.4,
  pointSize = 2.0,
  labSize = 5.0,
  labCol = 'black',
  labFace = 'plain',
  labhjust = 0.5,
  labvjust = 1.5,
  boxedLabels = FALSE,
  shape = 19,
  shapeCustom = NULL,
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  colCustom = NULL,
  colAlpha = 1/2,
  colGradient = NULL,
  colGradientBreaks = c(pCutoff, 1.0),
  colGradientLabels = c('0', '1.0'),
  colGradientLimits = c(0, 1.0),
  legendLabels = c('NS', expression(Log[2]~FC),
    'p-value', expression(p-value~and~log[2]~FC)),
  legendPosition = 'top',
  legendLabSize = 14,
  legendIconSize = 5.0,
  legendDropLevels = TRUE,
  encircle = NULL,
  encircleCol = 'black',
  encircleFill = 'pink',
  encircleAlpha = 3/4,
  encircleSize = 2.5,
  shade = NULL,
  shadeFill = 'grey',
  shadeAlpha = 1/2,
  shadeSize = 0.01,
  shadeBins = 2,
  drawConnectors = FALSE,
  widthConnectors = 0.5,
  typeConnectors = 'closed',
  endsConnectors = 'first',
  lengthConnectors = unit(0.01, 'npc'),
  colConnectors = 'grey10',
  arrowheads = TRUE,
  hline = NULL,
  hlineType = 'longdash',
  hlineCol = 'black',
  hlineWidth = 0.4,
  vline = NULL,
  vlineType = 'longdash',
  vlineCol = 'black',
  vlineWidth = 0.4,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  border = 'partial',
  borderWidth = 0.8,
  borderColour = 'black', 
  raster = FALSE)
{
  if(!is.numeric(toptable[[x]])) {
    stop(paste(x, ' is not numeric!', sep=''))
  }

  if(!is.numeric(toptable[[y]])) {
    stop(paste(y, ' is not numeric!', sep=''))
  }
  
  if (raster) {
    geom_point <- geom_point_rast
  }

  i <- xvals <- yvals <- Sig <- NULL

  toptable <- as.data.frame(toptable)
  toptable$Sig <- 'NS'
  toptable$Sig[(abs(toptable[[x]]) > FCcutoff)] <- 'FC'
  toptable$Sig[(toptable[[y]] < pCutoff)] <- 'P'
  toptable$Sig[(toptable[[y]] < pCutoff) &
    (abs(toptable[[x]]) > FCcutoff)] <- 'FC_P'
  toptable$Sig <- factor(toptable$Sig,
    levels=c('NS','FC','P','FC_P'))

  # some software programs return 0 for very low p-values
  # These throw an error in EnhancedVolcano
  # Detect these, issue warning, and convert these to
  # machine-lowest value possible
  #####
  # New functionality in > v1.2:
  # Now convert to 10^-1 lower than lowest non-zero p-value
  if (min(toptable[[y]], na.rm=TRUE) == 0) {
    # <= v1.2
    #warning(paste("One or more P values is 0.",
    #  "Converting to minimum possible value..."),
    #  call. = FALSE)
    #toptable[which(toptable[[y]] == 0), y] <- .Machine$double.xmin
    warning(paste('One or more p-values is 0.',
      'Converting to 10^-1 * current',
      'lowest non-zero p-value...'),
      call. = FALSE)
    toptable[which(toptable[[y]] == 0), y] <- min(
      toptable[which(toptable[[y]] != 0), y],
      na.rm = TRUE) * 10^-1
  }

  toptable$lab <- lab
  toptable$xvals <- toptable[[x]]
  toptable$yvals <- toptable[[y]]

  # If user has supplied values in selectLab, convert labels to
  # NA and then re-set with those in selectLab
  if (!is.null(selectLab)) {
    names.new <- rep(NA, length(toptable$lab))
    indices <- which(toptable$lab %in% selectLab)
    names.new[indices] <- toptable$lab[indices]
    toptable$lab <- names.new
  }

  # create a base theme that will later be modified
  th <- theme_bw(base_size = 24) +

    theme(
      legend.background = element_rect(),

      # title, subtitle, and caption
      plot.title = element_text(
        angle = 0,
        size = titleLabSize,
        face = 'bold',
        vjust = 1),
      plot.subtitle = element_text(
        angle = 0,
        size = subtitleLabSize,
        face = 'plain',
        vjust = 1),
      plot.caption = element_text(
        angle = 0,
        size = captionLabSize,
        face = 'plain',
        vjust = 1),

      # axis text
      axis.text.x = element_text(
        angle = 0,
        size = axisLabSize,
        vjust = 1),
      axis.text.y = element_text(
        angle = 0,
        size = axisLabSize,
        vjust = 0.5),
      axis.title = element_text(
        size = axisLabSize),

      # legend
      legend.position = legendPosition,
      legend.key = element_blank(),
      legend.key.size = unit(0.5, 'cm'),
      legend.text = element_text(
        size = legendLabSize),
      title = element_text(
        size = legendLabSize),
      legend.title = element_blank())

  # Create the plot object differently based on whether colCustom 
  # and shapeCustom are NULL or not. This helps to avoid messing up
  # the legend.
  #
  # 1, both colCustom and shapeCustom are activated
  if (!is.null(colCustom) & !is.null(shapeCustom)) {

    plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

      # over-ride legend icon sizes for colour and shape.
      # guide_legends are separate for colour and shape;
      # so, legends will be drawn separate
      guides(
        colour = guide_legend(
          order = 1,
          override.aes = list(
            size = legendIconSize)),
        shape = guide_legend(
          order = 2,
          override.aes = list(
            size = legendIconSize))) +

      # include new shape and colour encodings as aes
      geom_point(
        aes(
          color = factor(names(colCustom)),
          shape = factor(names(shapeCustom))),
        alpha = colAlpha,
        size = pointSize,
        na.rm = TRUE) +

      # specify the colour and shape with the supplied encoding
      scale_color_manual(values = colCustom) +
      scale_shape_manual(values = shapeCustom)

  # 2, only colCustom is activated and 'shape' has just a single value
  } else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 1) {

    plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

      # over-ride legend icon sizes for colour and shape.
      # guide_legends are separate for colour and shape;
      # so, legends will be drawn separate IF shape is also
      # included as aes to geom_point (it is not, here)
      guides(
        colour = guide_legend(
          order = 1,
          override.aes = list(
            size = legendIconSize)),
        shape = guide_legend(
          order = 2,
          override.aes = list(
            size = legendIconSize))) +

      # include new colour encodings as aes.
      # 'shape' is included, but outside aes
      geom_point(
        aes(
          color = factor(names(colCustom))),
        alpha = colAlpha,
        shape = shape,
        size = pointSize,
        na.rm = TRUE) +

      # specify the colour with the supplied encoding
      scale_color_manual(values = colCustom) +

      # 'shape' is not included as aes. Specifying guide = TRUE
      # here will result in legends merging
      scale_shape_manual(guide = TRUE)

  # 3, only colCustom is activated and 'shape' has 4 values
  } else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 4) {

    plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

      # over-ride legend icon sizes for colour and shape.
      # guide_legends are separate for colour and shape;
      # so, legends will be drawn separate
      guides(
        colour = guide_legend(
          order = 1,
          override.aes = list(
            size = legendIconSize)),
        shape = guide_legend(
          order = 2,
          override.aes = list(
            size = legendIconSize))) +

      # include new colour encodings as aes.
      # 'shape' is included in aes and mapped to 4
      # categories of NS, FC, P, FC_P
      geom_point(
        aes(
          color = factor(names(colCustom)),
          shape = Sig),
        alpha = colAlpha,
        size = pointSize,
        na.rm = TRUE) +

      # specify the colour with the supplied encoding
      scale_color_manual(values = colCustom) +

      # as it is included as aes, a separate legend
      # for 'shape' will be drawn. Here, over-ride that
      # legend
      scale_shape_manual(
        values = c(
          NS = shape[1],
          FC = shape[2],
          P = shape[3],
          FC_P = shape[4]),
        labels = c(
          NS = legendLabels[1],
          FC = legendLabels[2],
          P = legendLabels[3],
          FC_P = legendLabels[4]),
        guide = TRUE,
        drop = legendDropLevels)

  # 4, only shapeCustom is activated
  } else if (is.null(colCustom) & !is.null(shapeCustom)) {

    if (is.null(colGradient)) {

      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        # over-ride legend icon sizes for colour and shape.
        # guide_legends are separate for colour and shape;
        # so, legends will be drawn separate
        guides(
          colour = guide_legend(
            order = 1,
            override.aes = list(
              size = legendIconSize)),
          shape = guide_legend(
            order = 2,
            override.aes = list(
              size = legendIconSize))) +

        # include new shape encodings as aes.
        # Standard colour for NS, FC, P, FC_P,
        # are added to aes, too.
        geom_point(
          aes(
            color = Sig,
            shape = factor(names(shapeCustom))),
          alpha = colAlpha,
          size = pointSize,
          na.rm = TRUE) +

        # as it is included as aes, a separate legend
        # for 'colour' will be drawn. Here, over-ride that
        # legend
        scale_color_manual(
          values = c(
            NS = col[1],
            FC = col[2],
            P = col[3],
            FC_P = col[4]),
          labels = c(
            NS = legendLabels[1],
            FC = legendLabels[2],
            P = legendLabels[3],
            FC_P = legendLabels[4]),
          drop = legendDropLevels) +

        # specify the shape with the supplied encoding
        scale_shape_manual(values = shapeCustom)

    } else {

      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        # over-ride legend icon sizes for colour and shape.
        # guide_legends are separate for colour and shape;
        # so, legends will be drawn separate
        guides(
          shape = guide_legend(
            order = 2,
            override.aes = list(
              size = legendIconSize))) +

        # include new shape encodings as aes.
        # Standard colour for NS, FC, P, FC_P,
        # are added to aes, too.
        geom_point(
          aes(
            color = Sig,
            shape = factor(names(shapeCustom))),
          alpha = colAlpha,
          size = pointSize,
          na.rm = TRUE) +

        scale_colour_gradient(
          low = colGradient[1],
          high = colGradient[2],
          limits = colGradientLimits,
          breaks = colGradientBreaks,
          labels = colGradientLabels)

        # specify the shape with the supplied encoding
        scale_shape_manual(values = shapeCustom)

    }

  # 5, both colCustom and shapeCustom are null;
  # only a single shape value specified
  } else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 1) {

    if (is.null(colGradient)) {

      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        # over-ride legend icon sizes for colour and shape.
        # including 'shape' in the colour guide_legend here
        # results in the legends merging
        guides(colour = guide_legend(
          order = 1,
          override.aes = list(
            shape = shape,
            size = legendIconSize))) +

        geom_point(
          aes(color = Sig),
          alpha = colAlpha,
          shape = shape,
          size = pointSize,
          na.rm = TRUE) +

        scale_color_manual(
          values = c(
            NS = col[1],
            FC = col[2],
            P = col[3],
            FC_P = col[4]),
          labels = c(
            NS = legendLabels[1],
            FC = legendLabels[2],
            P = legendLabels[3],
            FC_P = legendLabels[4]),
          drop = legendDropLevels)

    } else {

      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        geom_point(
          aes(color = yvals),
          alpha = colAlpha,
          shape = shape,
          size = pointSize,
          na.rm = TRUE) +

        scale_colour_gradient(
          low = colGradient[1],
          high = colGradient[2],
          limits = colGradientLimits,
          breaks = colGradientBreaks,
          labels = colGradientLabels)
    }

  # 6, both colCustom and shapeCustom are null;
  # four shape values are specified
  } else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 4) {

    if (is.null(colGradient)) {
      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        # over-ride legend icon sizes for colour and shape.
        # including 'shape' in the colour guide_legend here
        # results in the legends merging
        guides(colour = guide_legend(
          order = 1,
          override.aes = list(
            shape = c(
              NS = shape[1],
              FC = shape[2],
              P = shape[3],
              FC_P = shape[4]),
            size = legendIconSize))) +

        geom_point(
          aes(
            color = Sig,
            shape = Sig),
          alpha = colAlpha,
          size = pointSize,
          na.rm = TRUE) +

        scale_color_manual(
          values = c(
            NS = col[1],
            FC = col[2],
            P = col[3],
            FC_P = col[4]),
          labels = c(
            NS = legendLabels[1],
            FC = legendLabels[2],
            P = legendLabels[3],
            FC_P = legendLabels[4]),
          drop = legendDropLevels) +

        scale_shape_manual(
          values = c(
            NS = shape[1],
            FC = shape[2],
            P = shape[3],
            FC_P = shape[4]),
          guide = FALSE,
          drop = legendDropLevels)

    } else {

      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + th +

        geom_point(
          aes(
            color = yvals,
            shape = Sig),
          alpha = colAlpha,
          size = pointSize,
          na.rm = TRUE) +

        scale_colour_gradient(
          low = colGradient[1],
          high = colGradient[2],
          limits = colGradientLimits,
          breaks = colGradientBreaks,
          labels = colGradientLabels) +

        scale_shape_manual(
          values = c(
            NS = shape[1],
            FC = shape[2],
            P = shape[3],
            FC_P = shape[4]),
          guide = FALSE,
          drop = legendDropLevels)

    }
  }

  # add more elements to the plot
  plot <- plot +

    xlab(xlab) +
    ylab(ylab) +

    xlim(xlim[1], xlim[2]) +
    ylim(ylim[1], ylim[2]) +

    geom_vline(xintercept = c(-FCcutoff, FCcutoff),
      linetype = cutoffLineType,
      colour = cutoffLineCol,
      size = cutoffLineWidth) +

    geom_hline(yintercept = -log10(pCutoff),
      linetype = cutoffLineType,
      colour = cutoffLineCol,
      size = cutoffLineWidth)

  # add elements to the plot for title, subtitle, caption
  plot <- plot + labs(title = title, 
    subtitle = subtitle, caption = caption)

  # add elements to the plot for vlines and hlines
  if (!is.null(vline)) {
    plot <- plot + geom_vline(xintercept = vline,
      linetype = vlineType,
      colour = vlineCol,
      size = vlineWidth)
  }
  if (!is.null(hline)) {
    plot <- plot + geom_hline(yintercept = -log10(hline),
      linetype = hlineType,
      colour = hlineCol,
      size = hlineWidth)
  }

  # Border around plot
  if (border == 'full') {
    plot <- plot + theme(panel.border = element_rect(
      colour = borderColour, fill = NA, size = borderWidth))
  } else if (border == 'partial') {
    plot <- plot + theme(axis.line = element_line(
      size = borderWidth, colour = borderColour),
      panel.border = element_blank(),
      panel.background = element_blank())
  } else {
    stop('Unrecognised value passed to \'border\'. Must be \'full\' or \'partial\'')
  }

  # Gridlines
  if (gridlines.major) {
    plot <- plot + theme(panel.grid.major = element_line())
  } else {
    plot <- plot + theme(panel.grid.major = element_blank())
  }
  if (gridlines.minor) {
    plot <- plot + theme(panel.grid.minor = element_line())
  } else {
    plot <- plot + theme(panel.grid.minor = element_blank())
  }

  # user has specified to draw with geom_text or geom_label?
  if (!boxedLabels) {

    # For labeling with geom_[text|label]_repel and
    # geom_[text|label] with check_overlap = TRUE, 4 possible
    # scenarios can arise
    if (drawConnectors && is.null(selectLab)) {

      if (arrowheads) {
        arr <- arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors)
      } else {
        arr <- NULL
      }

      plot <- plot + geom_text_repel(
        data = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff),
        aes(label = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff)[["lab"]]),
        size = labSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arr,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    } else if (drawConnectors && !is.null(selectLab)) {

      if (arrowheads) {
        arr <- arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors)
      } else {
        arr <- NULL
      }

      plot <- plot + geom_text_repel(
        data = subset(toptable,
          !is.na(toptable[['lab']])),
        aes(label = subset(toptable,
          !is.na(toptable[['lab']]))[['lab']]),
        size = labSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arr,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    } else if (!drawConnectors && !is.null(selectLab)) {

      plot <- plot + geom_text(
        data = subset(toptable,
          !is.na(toptable[['lab']])),
        aes(
          label = subset(toptable,
            !is.na(toptable[['lab']]))[['lab']]),
        size = labSize,
        check_overlap = TRUE,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    } else if (!drawConnectors && is.null(selectLab)) {

      plot <- plot + geom_text(
        data = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff),
        aes(label = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff)[['lab']]),
        size = labSize,
        check_overlap = TRUE,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)
    }

  } else {

    # For labeling with geom_[text|label]_repel and
    # geom_[text|label] with check_overlap = TRUE, 4 possible
    # scenarios can arise
    if (drawConnectors && is.null(selectLab)) {

      if (arrowheads) {
        arr <- arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors)
      } else {
        arr <- NULL
      }

      plot <- plot + geom_label_repel(
        data = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff),
        aes(label = subset(toptable,
          toptable[[y]]<pCutoff &
            abs(toptable[[x]]) > FCcutoff)[['lab']]),
        size = labSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arr,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    } else if (drawConnectors && !is.null(selectLab)) {

      if (arrowheads) {
        arr <- arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors)
      } else {
        arr <- NULL
      }

      plot <- plot + geom_label_repel(
        data = subset(toptable,
          !is.na(toptable[['lab']])),
        aes(label = subset(toptable,
          !is.na(toptable[['lab']]))[['lab']]),
        size = labSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arr,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    } else if (!drawConnectors && !is.null(selectLab)) {

      plot <- plot + geom_label(
        data = subset(toptable,
          !is.na(toptable[["lab"]])),
        aes(
          label = subset(toptable,
            !is.na(toptable[['lab']]))[['lab']]),
        size = labSize,
        #check_overlap = TRUE,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    } else if (!drawConnectors && is.null(selectLab)) {

      plot <- plot + geom_label(
        data = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff),
        aes(label = subset(toptable,
          toptable[[y]] < pCutoff &
            abs(toptable[[x]]) > FCcutoff)[['lab']]),
        size = labSize,
        #check_overlap = TRUE,
        hjust = labhjust,
        vjust = labvjust,
        colour = labCol,
        fontface = labFace,
        na.rm = TRUE)

    }
  }

  # encircle
  if (!is.null(encircle)) {
    plot <- plot + 
      geom_encircle(
        data = subset(toptable,
          rownames(toptable) %in% encircle),
        colour = encircleCol,
        fill = encircleFill,
        alpha = encircleAlpha,
        size = encircleSize,
        show.legend = FALSE,
        na.rm = TRUE)
  }

  # shade
  if (!is.null(shade)) {
    plot <- plot + 
      stat_density2d(
        data = subset(toptable,
          rownames(toptable) %in% shade),
        fill = shadeFill,
        alpha = shadeAlpha,
        geom = 'polygon',
        contour = TRUE,
        size = shadeSize,
        bins = shadeBins,
        show.legend = FALSE,
        na.rm = TRUE)
  }

  return(plot)
}

