# calculates upper and lower percentiles of lower depth
# adapted from https://github.com/KHanghoj/leopardpaper/blob/master/mapping_and_qc/reference_filter/depth/04_quantile_thresholds.R

sink(file(snakemake@log[[1]], open="wt"), type = "message")

library(ggplot2)

chrom.dp.table <- read.table(snakemake@input[[1]], header = FALSE)
genome.dp <- colSums(chrom.dp.table)
df <- data.frame(matrix(nrow = length(genome.dp), ncol = 0))
df$dp <- seq(from = 0, to = length(genome.dp)-1)
df$count <- genome.dp
mean <- with(df, mean(rep(x = dp, times = count)))
median <- with(df, median(rep(x = dp, times = count)))
genome.dp.cumsum <- cumsum(as.numeric(genome.dp))

if (snakemake@params[["method"]] == "percentile") {
  qup <- genome.dp.cumsum[length(genome.dp.cumsum)]*snakemake@params[["upper"]]
  qlow <- genome.dp.cumsum[length(genome.dp.cumsum)]*snakemake@params[["lower"]]
  upper <- min(which(genome.dp.cumsum > qup))
  lower <- min(which(genome.dp.cumsum > qlow))
  quants <- c(lower, upper)
} else if (snakemake@params[["method"]] == "median") {
  quants <- c(
    median * snakemake@params[["lower"]],
    median * snakemake@params[["upper"]]
  )
}

svg(file=snakemake@output[["plot"]])
xlim <- min(which(genome.dp.cumsum > genome.dp.cumsum[length(genome.dp.cumsum)]*0.995))
cov <- with(df, rep(x = dp, times = count))
hist(cov[cov <= xlim],
  main = paste0(
    "Global depth distribution ",
    snakemake@wildcards[["population"]],
    " depth samples"
  ),
  xlab = "Global Depth",
  ylab = "Count",
  breaks = 100
)
abline(v=quants[1], col = "red")
abline(v=quants[2], col = "red")
dev.off()

toprint <- c(mean,as.integer(quants),median)
write(toprint, snakemake@output[["summ"]])