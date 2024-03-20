print.usage <- function() {
	cat('\nUsage: plot_distance_count_multi.log.R <directory> [<directory> ...] <output>\n',file=stderr())
	cat('   <directory>, data directory (e.g., "CustardPyResults_Hi-C/Juicer_hg38/Sample")\n',file=stderr())
	cat('   <output>, output filename (e.g., "plot.pdf") \n',file=stderr())
	cat('\n',file=stderr())
}

args <- commandArgs(trailingOnly=TRUE)
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  print.usage()
  q()
}
folders <- args[-length(args)]
outputfile <- tail(args, n = 1)

library(reshape2)
library(ggplot2)
library(dplyr)

combined_data <- data.frame()
for (folder in folders) {
  filepath <- file.path(folder, "/distance/distance_vs_count.MAPQ30.log.txt")

  if (file.exists(filepath) && file.info(filepath)$size > 0) {
    data <- read.table(filepath)
    data <- data[c(3,5)]
    names(data) <- c("distance", "count")
    data$Sample <- basename(folder)
    combined_data <- rbind(combined_data, data)
  }
}

combined_data <- combined_data %>%
  group_by(Sample) %>%
  mutate(total_cis_contact = sum(count),
         probability = count / total_cis_contact)

xbreaks <- c(10000, 100000, 1000000, 10000000, 100000000, 1000000000)

filtered_data <- combined_data[combined_data$distance > 10000, ]
max_y_value <- max(filtered_data$probability)

pdf(outputfile, width = 6, height = 6)
g <- ggplot(combined_data, aes(distance, probability, colour = Sample))
g <- g + scale_x_log10(breaks=xbreaks, labels=xbreaks, limits=c(10000, NA))
g <- g + scale_y_continuous(limits = c(NA, max_y_value))
g <- g + geom_line()
g <- g + theme_minimal() +
  theme(panel.background = element_rect(fill = "white", colour = NA))
print(g)
dev.off()
