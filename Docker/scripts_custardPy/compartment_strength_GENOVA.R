library(GENOVA)
library(data.table) 
library(tidyr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
coolfile <- args[1]
sample_name <- args[2]
#resolution <- as.numeric(args[3])
#Peak <- args[3]

loaded_contacts <- load_contacts(
    signal_path = coolfile,
    sample_name = sample_name,
    balancing = TRUE,
    colour = "black" 
)

#H3K27acPeaks = read.delim(Peak, header = FALSE)
#CS_out = compartment_score(loaded_contacts, bed = H3K27acPeaks)
CS_out = compartment_score(loaded_contacts)

saddle_out = saddle(loaded_contacts, CS_discovery = CS_out, bins = 50)
CSS <- quantify(saddle_out)
compared <- tidyr::spread(unique(CSS[,-c(3,4)]), key = 'exp', value = 'strength')
compared_filtered <- compared %>% filter(chr != "chrYp")
col_means <- colMeans(compared_filtered[,-1], na.rm = TRUE)

write.table(col_means, 
            file = paste0("GENOVA_compartment_score.", sample_name, ".txt"), 
            quote = FALSE, row.names = TRUE, col.names = FALSE, sep = "\t")