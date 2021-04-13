suppressMessages(library(gridExtra))
suppressMessages(library(grid))
suppressMessages(library(ggplot2))
suppressMessages(library(lattice))
suppressMessages(library(stringr))
suppressMessages(library(ggpubr))

df <-read.csv("../results/corrplot_d23vd22.csv", header = TRUE, check.names = FALSE) 
df <- df %>% select(-patient_id)

colnames(df) <- gsub("-","", colnames(df))

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

imageWidth = 1300
imageHeight = 1300
dpi = 300

png(
  filename="../plots/scatter_plot_v2.png",
  width=imageWidth,
  height=imageHeight,
  units="px",
  pointsize=4,
  bg="white",
  res=dpi,
  type="cairo")

p <- grid.arrange(grobs = grob.list[select_grobs(lay)], layout_matrix = lay)

p

dev.off()
