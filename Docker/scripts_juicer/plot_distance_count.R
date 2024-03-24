print.usage <- function() {
	cat('\nUsage: plot_distance_count.R <input> <output> <label>\n',file=stderr())
	cat('   <input>, input data (e.g., "<CustardPyDir>/distance/distance_vs_count.10kb.MAPQ30.txt")\n',file=stderr())
	cat('   <output>, output filename (e.g., "plot.pdf") \n',file=stderr())
	cat('   <label>, sample label (e.g., "Sample1") \n',file=stderr())
	cat('\n',file=stderr())
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  print.usage()
  q()
}
inputfile <- args[1]
outputfile <- args[2]
label <- args[3]

library(reshape2)
library(ggplot2)
input1 <- read.table(inputfile)
data1 <- input1[c(3,5)]

total_cis_contact1 <- apply(data1[c(2)], 2, sum)
prob1 <- sapply(data1[c(2)], function(x){return(x/total_cis_contact1)})
p1 <- as.data.frame(prob1)
problist <- cbind(data1[c(1)],p1[c(1)])
name <- c("distance",label)
names(problist) <- name
df <- melt(problist, id.vars = "distance", measure.vars = c(label), variable.name = "data", value.name = "probability")

xbreaks <- c(10000, 100000, 1000000, 10000000, 100000000, 1000000000)
ybreaks <- c(1, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001)

g <- ggplot(df, aes(distance, probability, colour = data))
g <- g + scale_x_log10(breaks=xbreaks, labels=xbreaks)
g <- g + scale_y_log10(breaks=ybreaks, labels=ybreaks)
g <- g + geom_line()

# Modifying the theme for a plain background
g <- g + theme_minimal() +
  theme(panel.background = element_rect(fill = "white", colour = NA))

# Plotting
pdf(outputfile, width = 6, height = 5)
print(g)
dev.off()
