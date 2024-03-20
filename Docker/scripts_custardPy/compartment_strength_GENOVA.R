print.usage <- function() {
	cat('\nUsage: compartment_strength_GENOVA.R <coolfile> <sample_name>\n',file=stderr())
	cat('\n',file=stderr())
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  print.usage()
  q()
}
coolfile <- args[1]
sample_name <- args[2]
#resolution <- as.numeric(args[3])
#Peak <- args[3]

coolfile
sample_name
#Peak

library(GENOVA)
library(data.table)
library(tidyr)
library(dplyr)

loaded_contacts <- load_contacts(
    signal_path = coolfile,
    sample_name = sample_name,
    balancing = TRUE,
    colour = "black"
)

#peaks <- read.delim(Peak, header = FALSE)
#CS_out <- compartment_score(loaded_contacts, bed = peaks)
CS_out = compartment_score(loaded_contacts)

saddle_out <- saddle(loaded_contacts, CS_discovery = CS_out, bins = 50)
CSS <- quantify(saddle_out)
compared <- tidyr::spread(unique(CSS[,-c(3,4)]), key = 'exp', value = 'strength')
#print(compared)
compared_filtered <- compared %>% filter(chr != "chrYp")
col_means <- colMeans(compared_filtered[2], na.rm = TRUE)
#print(col_means)
write.table(col_means,
            file = paste0(sample_name, ".GENOVA_compartment_score.txt"),
            quote = FALSE, row.names = TRUE, col.names = FALSE, sep = "\t")