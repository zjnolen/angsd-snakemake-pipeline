sink(file(snakemake@log[[1]], open="wt"), type = "message")

aggregate_heterozygosity <- function(sfslist, pop) {

  s <- as.data.frame(read.table(pop, header = TRUE))

  hz <- data.frame(matrix(NA, ncol = 3, nrow = 0))
  colnames(hz) <- c("sample","pop","heterozygous_per_1000bp")

  for (i in 1:length(sfslist)) {
    
    sfs <- scan(sfslist[i])
    samplehz <- sfs[2]/sum(sfs)
    hz[nrow(hz) + 1,] <- c(s$sample[i],s$population[i],samplehz*1000)
  }

  hz$pop <- as.factor(hz$pop)
  hz$pop <- factor(hz$pop, levels = sort(unique(hz$pop)))
  hz[,3] <- as.numeric(hz[,3])
  return(hz)

}

bootstrap_heterozygosity <- function (bootlist, pop) {

  s <- as.data.frame(read.table(pop, header = TRUE))

  hz <- data.frame(matrix(NA, ncol = 3, nrow = 0))
  colnames(hz) <- c("sample", "lower_95_CI", "upper_95_CI")

  for (i in 1:length(bootlist)) {
    
    sfs <- read.table(bootlist[i], header = FALSE)
    sampleboothz <- sfs$V2/(sfs$V1+sfs$V2+sfs$V3)
    ci <- quantile(sampleboothz, probs = c(0.025, 0.975), names = FALSE)*1000
    hz[nrow(hz) + 1,] <- c(s$sample[i],ci[1],ci[2])
  }

  hz[,2] <- as.numeric(hz[,2])
  hz[,3] <- as.numeric(hz[,3])
  return(hz)
}

plot_heterozygosity <- function(agghz, popplot, indplot) {
  
  require(ggplot2)

  ggplot(agghz, aes(x=pop, y=heterozygous_per_1000bp, fill=pop)) +
    geom_boxplot() +
    xlab("Population") +
    ylab("Heterozygous sites per 1000bp") +
    theme_classic()
  
  ggsave(popplot)
  
  ggplot(agghz, aes(x=sample, y=heterozygous_per_1000bp)) +
    xlab("Sample") +
    ylab("Heterozygous sites per 1000bp") +
    geom_point() +
    facet_grid(cols = vars(pop), scales = "free_x") +
    geom_errorbar(aes(x=sample, ymin=lower_95_CI, ymax=upper_95_CI)) +
    theme_classic()

  ggsave(indplot)

}

agghz <- aggregate_heterozygosity(
  snakemake@input[["sfs"]],snakemake@input[["popfile"]])

boothz <- bootstrap_heterozygosity(
  snakemake@input[["bootsfs"]],snakemake@input[["popfile"]]
)

allhet <- merge(agghz, boothz, by = "sample")

write.table(allhet, file = snakemake@output[["table"]], row.names=FALSE,
  sep = "\t", quote=FALSE)

plot_heterozygosity(allhet, snakemake@output[["popplot"]],
  snakemake@output[["indplot"]])