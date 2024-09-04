## ------------------------------------------------------------------------
#This script adds taxonomy to each of the inferred ASVs

library(dada2)
library(Biostrings)
library(tidyverse)

cat(paste0('\n',"You are using DADA2 version ", packageVersion('dada2'),'\n'))

cat('################################\n\n')

args <- commandArgs(trailingOnly = TRUE)

seqtab.nochim <- args[1]
output <- args[2]
name <- args[3]
tax_db <- strsplit(args[4], ",")[[1]]
method <- args[5]
threshold <- as.integer(args[6])

dir.create(file.path(output, "03_taxonomy"), showWarnings = FALSE)
dir.create(file.path(output, "03_taxonomy", name), showWarnings = FALSE)

output <- paste0(output,"/03_taxonomy/",name,"/")

# Assign taxonomy (general)
set.seed(42) #random  generator necessary for reproducibility

seqtab.nochim <- readRDS(seqtab.nochim)

if (grepl('[Dd]ecipher|DECIPHER', method)){ # use decipher

  library(DECIPHER)
  cat('# The taxonomic classification included in decipher will be used\n')
  cat(paste0("# You are using DECIPHER version ", packageVersion('DECIPHER'),'\n'))
  cat('################################\n')

  dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs

  if (grepl('\\.RData$', tax_db[1])){
    load(tax_db[1])
  } else if (grepl('\\.rds$',tax_db[1])){
    trainingSet <- readRDS(tax_db[1])
  }

  if (is.na(threshold)){
      threshold <- 60
  }

  ids <- IdTaxa(dna,
              trainingSet,
              strand = "top",
              processors = NULL,
              verbose = FALSE,
              threshold = threshold)

  ranks <- paste0('level',1:max(trainingSet$levels))

  # Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy

  taxid <-
    names(ids) %>%
    map_df(~ tibble(seq_name = .x, tax = ids[[.x]]$taxon, ranks = ranks[1:length(ids[[.x]]$taxon)])) %>%
    pivot_wider(names_from = ranks, values_from=tax)

  write_tsv(taxid, paste0(output, name, "_tax_assignation_decipher.txt"))

  # Convert learnTaxa to data frame too

  n_seq <- length(ids)
  df_rows <- list()

  for(i in 1:n_seq){
    seq_name <- names(ids[i])
    taxonomy<- ids[[i]]$taxon
    confidence <- ids[[i]]$confidence
    df_rows[[i]] = data.frame(seq_name, taxonomy, confidence, taxo_level=ranks[1:length(ids[[i]]$taxon)])
  }

  df <- purrr::reduce(df_rows, bind_rows) %>%
    filter(taxo_level != "Root") %>%
    pivot_wider(names_from = taxo_level, values_from = c(taxonomy, confidence))

    write_tsv(df, paste0(output, name, "_decipher_output.txt"))

} else { # use regular dada2 classificator
  
  if (is.na(threshold)){
   threshold <- 50
  }

  cat('# The taxonomic classification included in dada2 will be used\n')
  cat(paste0('\n',"# You are using DADA2 version ", packageVersion('dada2'),'\n'))
  cat('################################\n')

  n_ranks <- length(str_split(str_remove(names(readDNAStringSet(tax_db[1], nrec = 1)),';$'),';')[[1]])
  ranks <- paste0('level',1:n_ranks)

  taxid <- assignTaxonomy(seqtab.nochim,
                        tax_db[1], 
                        multithread=TRUE,
                        minBoot=threshold,
                        taxLevels = ranks,
                        outputBootstraps = TRUE)
  cat('# Taxonomy assigned to genus level\n')
  
  if (!is.na(tax_db[2])) { # add species level if db available
    taxid <- addSpecies(taxid$tax, 
                        tax_db[2], 
                        verbose=TRUE, 
                        allowMultiple=3)
    cat('\n# Taxonomy assigned to species level\n')    
  }

  rownames(taxid$tax) <- colnames(seqtab.nochim)
  rownames(taxid$boot) <- colnames(seqtab.nochim)

  tax_table <-
    taxid$tax %>%
    as_tibble(rownames = 'seq_name')
  
  boot_table <-
    taxid$boot %>%
    as_tibble(rownames = 'seq_name')
  
  write_tsv(tax_table, paste0(output, name, "_tax_assignation_dada.txt"))
  write_tsv(boot_table, paste0(output, name, "_dada_output.txt"))

}

# Farewell
cat('\n')
cat(paste0('# The obtained taxonomy file can be found in "', paste0(output, name, "_tax_assignation_",method,".txt"), '"\n'))
cat('\n# All done!\n\n')
