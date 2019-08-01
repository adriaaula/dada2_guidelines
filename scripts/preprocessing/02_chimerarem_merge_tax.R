## ------------------------------------------------------------------------
#This script merges together all tables generated in the runs and assigns taxonomy to
#each ASV from the whole table

library(dada2)
library(tidyverse)

cat(paste0('\n',"You are using DADA2 version ", packageVersion('dada2'),'\n\n'))

cat('################################\n\n')

args <- commandArgs(trailingOnly = TRUE)

seqtables <- strsplit(args[1], ",")[[1]]
output <- args[2]
name <- args[3]
tax_db <- args[4]
tax_db_sp <- args[5]
trim_length <- as.integer(strsplit(args[6], ",")[[1]])

dir.create(file.path(output, "02_mergeruns-taxonomy"), showWarnings = FALSE)
dir.create(file.path(output, "02_mergeruns-taxonomy", name), showWarnings = FALSE)

output <- paste0(output,"/02_mergeruns-taxonomy/",name,"/")

# Get all the rds files
list.df <- map(seqtables, readRDS)

st.all <- mergeSequenceTables(tables = list.df)

track.final <- data.frame(sample = rownames(st.all),
                          raw = rowSums(st.all))

seqtab.raw <- removeBimeraDenovo(st.all, method="consensus",
                                     multithread=TRUE)

num.chimera.removed <- ncol(st.all) - ncol(seqtab.raw)
perc.num.chimera.removed <- round(100*num.chimera.removed/ncol(st.all) %>% round(2),2)
reads.chimera.removed <- sum(colSums(st.all)) - sum(colSums(seqtab.raw))
perc.reads.chimera.removed <- round(100*reads.chimera.removed/sum(colSums(st.all)) %>% round(2),2)

cat(paste0('# ',num.chimera.removed," chimera were found and removed\n"))
cat(paste0('# These represent ',perc.num.chimera.removed,'% of total ASVs and ',perc.reads.chimera.removed,'% of total reads\n\n'))

total <- sum(colSums(seqtab.raw))

track.final$chimera.removed  <- rowSums(seqtab.raw)

# Distribution of variants
cat("# The variants (ASVs) have the following length distribution:\n")
table(nchar(getSequences(seqtab.raw)))

# Trim the unespecific amplifications from our dataset
seqtab <- seqtab.raw[,nchar(colnames(seqtab.raw)) %in% seq(trim_length[1],
                                                           trim_length[2])]

cat(paste0('\n# Reads shorter than ',trim_length[1],'bp and longer than ',trim_length[2], 'bp were removed\n'))

# Afterwards:
cat("\n# The variants (ASVs) after length filtering have the following length distribution:\n")

table(nchar(getSequences(seqtab)))

track.final <- track.final %>%
                    mutate( too_long_variants = rowSums(seqtab))

final <- sum(colSums(seqtab))

cat(paste0("\n# A total of ", round((final * 100) / total, digits =2), "% reads were kept!\n\n"))

cat("# Collapsing at 100% id\n\n")

seqtab <- collapseNoMismatch(seqtab,  minOverlap = 50)

track.final <- track.final %>%
                    mutate( collapsed_100 = rowSums(seqtab))


saveRDS(seqtab, paste0(output, name, "_seqtab_final.rds"))

# Extract the FASTA from the dataset
uniquesToFasta(seqtab,
               paste0(output, name, "_seqtab_final.fasta"),
               ids= paste0("asv",c(1:ncol(seqtab)), ";size=", colSums(seqtab)))

# Assign taxonomy (general)
set.seed(42) #random  generator necessary for reproducibility

tax <- assignTaxonomy(seqtab,
                    tax_db,
                    multithread=TRUE,minBoot=80)
# minboot: N of idntical bootstraps to gen. identification

# Assign species
tax.sp <- addSpecies(tax, tax_db_sp, verbose=TRUE, allowMultiple=3)

# Write to disk
saveRDS(tax.sp, paste0(output, name, "_tax_assignation.rds"))


track.final <- track.final %>%
               mutate( diff.total = round(collapsed_100 / raw, digits =2))

write_tsv(track.final, paste0(output, name, "_track_analysis_final.tsv"))

cat(paste0('\n\n# Your final ASV table can be found in ', paste0(output, name, "_seqtab_final.rds"),'\n'))
cat(paste0('# A FASTA file with your final ASVs was written in ',paste0(output, name, "_seqtab_final.fasta"), '\n'))
cat(paste0('# The taxonomy file can be found in ', paste0(output, name, "_tax_assignation.rds"), '\n'))
cat('# In case you need it, we created a .txt file with tax and counts tables merged\n')
cat('\n# All done!\n\n')
