sink(file(snakemake@log[[1]], open="wt"), type = "message")

plot_admix <- function(qopt, pop, k, plotout) {
	q <- read.table(qopt)
	s <- as.data.frame(read.table(pop, header = TRUE, sep = "\t"))
	p <- s$population
	n <- nrow(s)

	ord <- order(p)

	svg(plotout, height = 5, width = 2+n*0.3)
	par(mar=c(8,4.1,4.1,2.1))
	barplot(t(q)[,ord],
		col=1:as.integer(k),
		names=p[ord],
		las=2,
		space=0,
		border=NA,
		ylab=paste0("Admixture proportions for K=",toString(k))
		)
	#mtext(text="Individuals", side=1, line=6.5)
	dev.off()
}

plot_admix(snakemake@input[[1]], snakemake@input[[2]],
	snakemake@wildcards[["kvalue"]], snakemake@output[[1]])