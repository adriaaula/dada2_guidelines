## ------------------------------------------------------------------------
#This script adds taxonomy to each of the inferred ASVs

library(dada2)
library(tidyverse)
library(DECIPHER)

cat(paste0('\n',"You are using DADA2 version ", packageVersion('dada2'),'\n'))
cat(paste0("You are using DECIPHER version ", packageVersion('DECIPHER'),'\n\n'))

cat('################################\n\n')

args <- commandArgs(trailingOnly = TRUE)

seqtab.nochim <- args[1]
output <- args[2]
name <- args[3]
tax_db <- args[4]
threshold <- as.integer(args[5])

dir.create(file.path(output, "03_taxonomy"), showWarnings = FALSE)
dir.create(file.path(output, "03_taxonomy", name), showWarnings = FALSE)

output <- paste0(output,"/03_taxonomy/",name,"/")

# Assign taxonomy (general)
set.seed(42) #random  generator necessary for reproducibility

seqtab.nochim <- readRDS(seqtab.nochim)

dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs

load(tax_db)

ids <- IdTaxa(dna,
             trainingSet,
             strand = "top",
             processors = NULL,
             verbose = FALSE,
             threshold = threshold)

ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest

# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
        m <- match(ranks, x$rank)
        taxa <- x$taxon[m]
        taxa[startsWith(taxa, "unclassified_")] <- NA
        taxa
}))

colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

# Write to disk
saveRDS(taxid, paste0(output, name, "_tax_assignation.rds"))

# Create a merged table with counts and tax

taxid <- as_tibble(taxid, rownames = 'ASV')
merged_table <- as_tibble(t(seqtab.nochim), rownames = 'ASV') %>%
                left_join(taxid, by = 'ASV') %>%
                mutate(ASV_id = paste0("asv",c(1:nrow(.)))) %>%
                select(ASV_id, everything())

write_tsv(x = merged_table,
          path = paste0(output, name, "_merged_table.txt"))

cat(paste0('# The obtained taxonomy file can be found in "', paste0(output, name, "_tax_assignation.rds"), '"\n'))
cat(paste0('# Although we always recommend you to work directly in R with .rds files, we created a .txt in "',paste0(output, name, "_merged_table.txt",'" with tax and counts tables merged\n')))
cat('\n# All done!\n\n')
