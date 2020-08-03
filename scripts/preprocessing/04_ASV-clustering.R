library(tidyverse)
library(DECIPHER)
library(Biostrings)

args <- commandArgs(trailingOnly = TRUE)

seqtab <- readRDS(args[1])
id <- (100 - as.integer(args[2]))/100
output <- args[3]
name.run <- args[4]

cat(output)
cat(name.run)
cat(id)
cat(args[2])
cat(class(args[2]))
cat(class(id))

asv_sequences <- colnames(seqtab)
sample_names <- rownames(seqtab)
dna <- Biostrings::DNAStringSet(asv_sequences)
asv_sizes <- colSums(seqtab)

## Find clusters of ASVs to form the new OTUs
alignment <- 
  DECIPHER::AlignSeqs(dna, 
                      processors = NULL)

dist_matrix <- 
  DECIPHER::DistanceMatrix(alignment,
                           processors = NULL)

clusters <- 
  DECIPHER::IdClusters(
  dist_matrix, 
  method = "complete",
  cutoff = id, # use `cutoff = 0.03` for a 97% OTU 
  processors = NULL)

cat(paste0('# Clustering at ',args[2],'% identity\n'))

cat(paste0('# ',
           length(unique(clusters$cluster)),
           ' OTUs created out of ',
           length(asv_sequences),
           ' ASVs\n'))

clusters <-
  tibble(ASV = asv_sequences,
         cluster = clusters$cluster) %>%
  mutate(size = asv_sizes) %>%
  group_by(cluster) %>%
  mutate(OTU = ASV[size == max(size)][1]) # avoid having duplicates when length(max(size)) > 1

summary_clustering <- 
  clusters %>% 
  group_by(cluster) %>% 
  tally() %>% 
  pull(n) %>% 
  summary()

cat('# Summary of cluster sizes:\n')
print(summary_clustering)

cat('\n# Creating the clustered seqtab...\n')

merged_seqtab <- 
  seqtab %>% 
  t() %>%
  rowsum(clusters$cluster) %>%
  as_tibble(rownames = 'cluster') %>%
  mutate(cluster = as.integer(cluster)) %>% 
  left_join(clusters %>% select(OTU, cluster) %>% unique(),
            by = 'cluster') %>% 
  select(-cluster) %>% 
  column_to_rownames('OTU') %>% 
  t()

saveRDS(merged_seqtab, paste0(output,name.run,"_seqtab_clust",args[2],".rds"))
write_tsv(clusters, paste0(output,name.run,"_ASV2OTU_clust",args[2],".tsv"))

cat(paste0('# An OTU clustered table was created, you can find it in "',
           paste0(output,name.run,"_seqtab_clust",args[2],".rds"),
           '"\n'))

cat(paste0('# An OTU to ASV correspondence table was created, you can find it in "',
           paste0(output,name.run,"_ASV2OTU_clust",args[2],".tsv"),
           '"\n'))

cat('# All done!\n')
