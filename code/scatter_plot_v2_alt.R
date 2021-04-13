suppressMessages(library(gridExtra))
suppressMessages(library(grid))
suppressMessages(library(ggplot2))
suppressMessages(library(lattice))
suppressMessages(library(stringr))
suppressMessages(library(ggpubr))

df <-read.csv("../results/corrplot_d23vd22.csv", header = TRUE, check.names = FALSE) 
df <- df %>% select(-patient_id)

colnames(df) <- gsub("-","", colnames(df))

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

imageWidth = (1300/4)*6
imageHeight = (1300/4)*3
dpi = 300

png(
  filename="../plots/scatter_plot_v2_alt.png",
  width=imageWidth,
  height=imageHeight,
  units="px",
  pointsize=4,
  bg="white",
  res=dpi,
  type="cairo")

grob.list <- lapply(seq_along(plots),plot_cor)
select_grobs <- function(hlay) {
  id <- unique(c(t(hlay))) 
  id[!is.na(id)]
} 
p <- grid.arrange(grobs = grob.list[select_grobs(lay)], layout_matrix = lay)
dev.off()
