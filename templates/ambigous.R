#!/apps/eb/2020b/skylake/software/R/4.2.1-foss-2022a/bin/Rscript

library(dplyr)

input <- read.table("${best_guess_folder}/hla/R1_bestguess_G.txt", header=TRUE)
input\$Allele <- if_else(input\$perfectG == 1, input\$Allele, "NA")
        
write.table(input\$Allele, "${dataset}_perfectG.tsv", col.names = T, row.names = F, quote = F, sep = "\\t")
