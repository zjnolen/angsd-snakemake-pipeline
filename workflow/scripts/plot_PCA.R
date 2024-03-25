sink(file(snakemake@log[[1]], open="wt"), type = "message")

library(ggplot2)

plot_pca <- function(cov, pop, xpc, ypc, plotout){

  C <- as.matrix(read.table(cov))
  pop <- read.table(pop, header = TRUE, sep = "\t")

  e <- eigen(C)

  pc <- data.frame(e$vectors[,as.integer(xpc):as.integer(ypc)])
  pc$pop <- pop$population
  pc$ind <- pop$sample
  loadings <- round((e$values/sum(e$values))*100, 2)

  ggplot(pc, aes(x=X1,y=X2,col=pop)) +
    geom_point(size=3) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    labs(x = paste0("PC",toString(xpc)," (",loadings[as.integer(xpc)],"%)"), 
      y = paste0("PC",toString(ypc)," (",loadings[as.integer(ypc)],"%)"))

  ggsave(plotout)

}

plot_pca(snakemake@input[[1]], snakemake@input[[2]], 
  snakemake@wildcards[["xpc"]], snakemake@wildcards[["ypc"]], 
  snakemake@output[[1]])