## ------------------------------------------------------------------------
library(stringr)
library(tidyverse)
library(dada2)

cat(paste0('\n',"You are using DADA2 version ", packageVersion('dada2'),'\n\n'))

cat('################################\n\n')

args <- commandArgs(trailingOnly = TRUE)

##Variables to change for EACH RUN
path<- args[1] #where the fastq files are located
output <-args[2] #name of analysis directory.
name.run <- args[3] #Run name. A directory for the output results will be created
trunclen <- as.integer(strsplit(args[4], ",")[[1]])
maxee <- as.integer(strsplit(args[5], ",")[[1]])

minover <- as.integer(args[6]) ##Default value is 15, not a good option go below 10.
pool_method <- args[7]

pool_method <- 
    if (is.na(pool_method)){
        pool_method <- FALSE
    } else if (grepl('pool',pool_method)){
        pool_method <- TRUE
    } else if (grepl('pseudo',pool_method)){
        pool_method <- 'pseudo'
    } else {
        pool_method <- FALSE
    }

## ------------------------------------------------------------------------

fnFs <- sort(list.files(path, pattern = 'R1.fastq'))
fnRs <- sort(list.files(path, pattern= 'R2.fastq'))

if (length(fnFs) == 0){
  stop("It could be that your data directory is empty.
       Otherwise filenames do not seem to follow the pattern <sample-name>_R1.fastq
       Please check these and run again")
}

sample.names <-
  map_chr(.x = fnFs,
          .f = ~ gsub("(.*)_R1\\.fastq.*", "\\1",.x))

if (length(unique(sample.names)) != length(fnFs)){
  stop("Filenames do not seem to follow the pattern <sample-name>_R1.fastq
       Please check these and run again")
}

fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

## ------------------------------------------------------------------------
#Check and create a dir for filtered fastqs
dir.create(file.path(path, paste0("filtered_",name.run)), showWarnings = FALSE)

filt_path <- file.path(path, paste0("filtered_",name.run)) # Filtered fasta in filtered_<name.run>/ subdir
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

## ------------------------------------------------------------------------

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=trunclen,
              maxEE=maxee,  rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)

cat('# Filtering and trimming done with the following parameters:\n')
cat(paste0('# Forward pair: trimming at ',trunclen[1],'bp and max expected error ',maxee[1],'\n'))
cat(paste0('# Reverse pair: trimming at ',trunclen[2],'bp and max expected error ',maxee[2],'\n\n'))

## ------------------------------------------------------------------------
#Check and create a dir for output results

dir.create(file.path(output, "01_errors-output"), showWarnings = FALSE)
dir.create(file.path(output, "01_errors-output", name.run), showWarnings = FALSE)

output <- str_c(output,"/01_errors-output/",name.run,"/")

## ------------------------------------------------------------------------
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

## ------------------------------------------------------------------------
err.plotf <- plotErrors(errF, nominalQ=TRUE)
ggsave(str_c(output,"errors_",name.run,"_fwd.pdf"),plot=err.plotf, width = 9, height = 8)
err.plotr <- plotErrors(errR, nominalQ=TRUE)
ggsave(str_c(output,"errors_",name.run,"_rev.pdf"),plot=err.plotr, width = 9, height = 8)

cat('\n# Errors learnt\n')

## ------------------------------------------------------------------------
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

cat('\n# Dereplication done\n\n')

## ------------------------------------------------------------------------
if (pool_method == TRUE){
    cat('\n# Performing dada inference with pooling\n\n')
    } else if (pool_method == FALSE) {
        cat('\n# Performing dada inference with no pooling\n\n')
    } else {
        cat(paste0('\n# Performing dada inference with pooling method ',pool_method,'\n\n'))
}


dadaFs <- dada(derepFs, err=errF, multithread=TRUE,pool=pool_method)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE,pool=pool_method)

cat('\n# DADA2 algorithm performed\n\n')

## ------------------------------------------------------------------------
mergers <- mergePairs(dadaFs, derepFs,
                      dadaRs, derepRs,
                      minOverlap = minover)

cat('# Pairs were merged\n\n')

## ------------------------------------------------------------------------
seqtab <- makeSequenceTable(mergers)

cat(paste0('# Number of samples: ',dim(seqtab)[1], '\n'))
cat(paste0('# Number of detected variants (ASVs): ',dim(seqtab)[2]))
cat("\n# The variants (ASVs) have the following length distribution:\n")
table(nchar(getSequences(seqtab)))

saveRDS(seqtab, str_c(output,name.run,"_seqtab.rds"))
## ------------------------------------------------------------------------

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN),
               sapply(mergers, getN),
               rowSums(seqtab))

colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled")
rownames(track) <- sample.names

track <- track %>%
         data.frame() %>%
         rownames_to_column( var = 'sample') %>%
         mutate(diff1 = round(filtered/input, 2),
                diff2 = round(denoised/filtered, 2),
                diff3 = round(merged/denoised, 2),
                diff4 = round(tabled/merged, 2),
                diff.total = round(tabled/input, 2))

cat("\n# The median of reads kept is the following:\n\n")

summary(track$diff.total)

write_tsv(data.frame(track),str_c(output,name.run,"_track_analysis.tsv"))

# The scripts generates an artifactural pdf in the main directory.
# We are going to remove it. It is a bad solution, but a solution after all!
if (file.exists('Rplots.pdf')) file.remove('Rplots.pdf')

cat(paste0('\n# An ASV table was created, you can find in "',str_c(output,name.run,"_seqtab.rds"),'"\n'))
cat('# Remember that this table still contains chimeras. You have now to run script "02_run-chimerarem_merge.sh" to remove them and add taxonomy\n')
cat(paste0('# In "',str_c(output,name.run,"_track_analysis.tsv"),"\" you will find a table where you can check the loss of reads in each step. Check it out to see if everything's correct!",'\n'))
cat('\n# All done!\n\n')

# END
