suppressMessages(library(corrplot)) # 0.84
suppressMessages(library(tidyverse)) # 1.3.0
suppressMessages(library(colorspace)) # 1.4.1

# From qdap v2.4.3
# https://cran.r-project.org/web/packages/qdap/index.html
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

df <- read.csv("../results/naive_plot_counts.csv", header = TRUE, check.names = FALSE)

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


imageWidth = 1300
imageHeight = 1300
dpi = 300

png(
  filename="../plots/corrplot_d2vd1.png",
  width=imageWidth,
  height=imageHeight,
  units="px",
  pointsize=4,
  bg="white",
  res=dpi,
  type="cairo")

p <- corrplot(M, method = "ellipse", col = col,
              type = "upper", order = "hclust", number.cex = 1, number.font=2,
              addCoef.col = "black", # Add coefficient of correlation
              tl.col = "black", tl.srt = 90, tl.cex=1.5, # Text label color and rotation
              # Combine with significance
              p.mat = p.mat, sig.level = 0.01, insig = "blank", 
              # hide correlation coefficient on the principal diagonal
              diag = FALSE,cl.cex = 2)

p

dev.off()

df <- df %>% rownames_to_column("patient_id")
write.csv(df, "../results/corrplot_d2vd1.csv", row.names = FALSE, quote = FALSE)