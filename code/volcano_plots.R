
  suppressMessages(library(stringr)) #1.4.0
  suppressMessages(library(tidyverse)) #1.3.0
  suppressMessages(library(EnhancedVolcano)) #1.6.0

  df <- read.csv("../results/naive_deg_analysis.csv", header = TRUE)
  
  contrasts <- c("d2_d1","d8_d1","d23_d22","d23_d22_d2_d1")
  subtitles <- c("Vaccination 1 ( 1 dpv vs pre)",
                "Vaccination 1 ( 7 dpv vs pre)",
                "Vaccination 2 ( 1 dpv vs pre)",
                "Vaccination 2 vs Vaccination 1 (1 dpv vs pre)")
  
  for(i in seq_along(contrasts)){
    c <- contrasts[i]
    label.col <- "Gene"
    sig.col <- paste0(c,"_adjpval")
    lfc.col <-paste0(c,"_logFC")
    columns_of_interest <- c(label.col,lfc.col,sig.col)
    
    df.sub <- df %>% dplyr::select(one_of(columns_of_interest)) %>% dplyr::filter(!is.na(!!rlang::sym(lfc.col))) #seurat introduces NAs
    
    change_lfc_name <- "log2FC"
    change_sig_name <- "FDR"
    colnames(df.sub) <- c(label.col,change_lfc_name,change_sig_name)
    
    cat(paste0("Genes in initial dataset: ", nrow(df.sub),"\n"))
    
    #Select top genes by logFC or Significance
    no_genes_to_label <- 30
    value_to_sort_the_output_dataset <- "p-value"
    if (value_to_sort_the_output_dataset=="fold-change") {
      df.sub %>% dplyr::arrange(desc(abs(!!rlang::sym(change_lfc_name)))) -> df.sub
    } else if (value_to_sort_the_output_dataset=="p-value") {
      df.sub %>% dplyr::arrange(!!rlang::sym(change_sig_name)) -> df.sub
    }
    genes_to_label <- as.character(df.sub[1:no_genes_to_label,label.col])
    
    additional_labels <- ""
    additional_labels <- unlist(str_split(additional_labels,","))
    filter <- additional_labels %in% df.sub[,label.col]
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
    subtitle <- subtitles[i]

    pCutoff  = 0.05
    FCcutoff = 1.0
    
    significant=as.vector(table(abs( df.sub[,change_lfc_name] ) > FCcutoff &
                                  df.sub[,change_sig_name]   < pCutoff))[2]
    
    # fix pvalue == 0
    shapeCustom <- rep(19,nrow(df.sub))
    maxy <-  max(-log10(df.sub[[change_sig_name]]), na.rm=TRUE)
    if(0 > 0){
      maxy <- 0
    }
    
    cat(paste0("Maxy: ",maxy,"\n"))
    if(maxy == Inf){
      # Sometimes, pvalues == 0
      keep <- df.sub[[change_sig_name]] > 0
      df.sub[[change_sig_name]][!keep] <- min(df.sub[[change_sig_name]][keep])
      shapeCustom[!keep] <- 17
      
      maxy <- -log10(min(df.sub[[change_sig_name]][keep]))
      cat("Some p-values equal zero. Adjusting y-limits.\n")
      cat(paste0("Maxy adjusted: ",maxy,"\n"))
      
    }
    
    # By default, nothing will be greater than maxy. User can set this value lower
    keep <- -log10(df.sub[[change_sig_name]]) <= maxy
    df.sub[[change_sig_name]][!keep] <- maxy
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
    
    filename = file.path("../plots/",paste0(c,"_volcano.png"))
    png(
      filename=filename,
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")
    
    p <- EnhancedVolcano(   df.sub,x=change_lfc_name,y=change_sig_name,
                            lab=df.sub[,label.col],
                            selectLab = genes_to_label,
                            title=paste0(title," (Significant=",significant,")"),
                            subtitle <- subtitle,
                            xlab=xlab,
                            ylab=bquote(~-Log[10]~.(change_sig_name)),
                            xlim=c(floor(min(df.sub[,change_lfc_name])) - xlim_additional,ceiling(max(df.sub[,change_lfc_name]))+ xlim_additional),
                            ylim=c(0, maxy + ylim_additional),
                            pCutoff=pCutoff,
                            FCcutoff=FCcutoff,
                            axisLabSize=axisLabSize,
                            labSize=labSize,
                            pointSize=pointSize,
                            shapeCustom=shapeCustom
    )
    print(p)
    dev.off()
  }
