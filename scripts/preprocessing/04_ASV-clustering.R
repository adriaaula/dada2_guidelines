suppressMessages(library(tidyverse))
suppressMessages(library(DECIPHER))
suppressMessages(library(Biostrings))
suppressMessages(library(data.table))

cat(paste0("# You are using DECIPHER version ", packageVersion('DECIPHER'),'\n'))

## read args

args <- commandArgs(trailingOnly = TRUE)

seqtab <- readRDS(args[1])
id <- as.numeric(args[2])
cutoff <- (100-id)/100
output <- args[3]
name.run <- args[4]
representative <- args[5]
min.coverage <- as.numeric(args[6])

dir.create(output, showWarnings = FALSE, recursive = T) # create output dir in case it does not exist

cat(paste0('# Clustering at ',id,'% identity\n'))

## Create ASV df

asv_df <- 
  tibble(
    seq = colnames(seqtab),
    size = colSums(seqtab)
  ) |> 
  mutate(length = nchar(seq),
         name = paste0('ASV_',row_number())) # add names to avoid duplicates

## Create DNAStringSet for clustering

dna <-
  asv_df |> 
  select(name, seq) |> 
  deframe() |> 
  DNAStringSet()

## Find clusters of ASVs 

clusters <-
  Clusterize(
    myXStringSet = dna,
    cutoff = cutoff,
    minCoverage = min.coverage,
    processors = NULL
  ) |> 
  as_tibble(rownames = 'name')

cat(paste0('# ',
           n_distinct(clusters$cluster),
           ' clusters created out of ',
           length(dna),
           ' ASVs\n'))

summary_clustering <- 
  clusters |> 
  group_by(cluster) |> 
  tally() |> 
  pull(n) |>  
  summary()

cat('# Summary of cluster sizes:\n')
print(summary_clustering)

asv_df_clusters <- 
  asv_df |> 
  left_join(clusters, by  = 'name')

if (representative == 'abundance'){
  
  cat('# Selecting cluster representatives by highest abundance\n')
  
  clusters_out <- 
    asv_df_clusters |> 
    arrange(-size, -length) |> 
    group_by(cluster) |> 
    mutate(representative = dplyr::first(seq)) |> 
    select(-name)
    
} else if (representative == 'length'){
  
  cat('# Selecting cluster representatives by highest length\n')
  
  clusters_out <- 
    asv_df_clusters |> 
    arrange(-length, -size) |> 
    group_by(cluster) |> 
    mutate(representative = dplyr::first(seq)) |> 
    select(-name)
  
} else {
  stop("Representative method selection must be one of the following: 'abundance', 'length'")
}

cat(paste0('# A correspondence table between ASVs and cluster representatives was created, you can find it in "',
           paste0(output,name.run,"_clusters",id,"id.tsv"),
           '"\n'))
write_tsv(clusters_out, paste0(output,name.run,"_clusters_",id,"id.tsv"))

cat('\n# Creating the clustered seqtab...\n')

representatives <- 
  clusters_out |> 
  ungroup() |> 
  select(seq, representative)

merged_seqtab <- 
  seqtab |> 
  as_tibble(rownames = 'sample') |> 
  pivot_longer(names_to = 'seq',
               values_to = 'abundance',
               -1) |> 
  left_join(representatives, by = 'seq') |> 
  group_by(representative, sample) |> 
  summarise(abundance = sum(abundance), .groups = 'drop_last') |> # add .groups to avoid annoying summarise warning
  pivot_wider(names_from = representative,
              values_from = abundance,
              values_fill = 0) |> 
  column_to_rownames('sample')
  
saveRDS(merged_seqtab, paste0(output,name.run,"_seqtab_clust_",id,"id.rds"))


cat(paste0('# A clustered table was created, you can find it in "',
           paste0(output,name.run,"_seqtab_clust",id,"id.rds"),
           '"\n'))



cat('# All done!\n')
