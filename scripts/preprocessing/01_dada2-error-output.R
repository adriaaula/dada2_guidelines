## ------------------------------------------------------------------------
library(stringr)
library(tidyverse)
library(dada2); packageVersion("dada2")

args <- commandArgs(trailingOnly = TRUE)

##Variables to change for EACH RUN
path<- args[1] #where the fastq files are located
output <-args[2] #name of analysis directory.
name.run <- args[3] #Run name. A directory for the output results will be created
trunclen <- as.integer(strsplit(args[4], ",")[[1]])
maxee <- as.integer(strsplit(args[5], ",")[[1]])

minover <- as.integer(args[6]) ##Default value is 15, not a good option go below 10. 

## ------------------------------------------------------------------------
head(list.files(path))

## ------------------------------------------------------------------------
fnFs <- sort(list.files(path, pattern="R1.fastq"))
fnRs <- sort(list.files(path, pattern="R2.fastq"))
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1) #maybe if there is no similar delim the script is not working here. Beware!

fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)


## ------------------------------------------------------------------------
#Check and create a dir for filtered fastqs
dir.create(file.path(path, "filtered"), showWarnings = FALSE)

filt_path <- file.path(path, "filtered") # Filtered fast in filtered/ subdir
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))

## ------------------------------------------------------------------------
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=trunclen,
              maxEE=maxee,  rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) 

head(out)

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

## ------------------------------------------------------------------------
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

## ------------------------------------------------------------------------
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)


## ------------------------------------------------------------------------
mergers <- mergePairs(dadaFs, derepFs,
                      dadaRs, derepRs,
                      minOverlap = minover)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

## ------------------------------------------------------------------------
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))


saveRDS(seqtab, str_c(output,name.run,"seqtab.rds"))
## ------------------------------------------------------------------------

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN),
               sapply(mergers, getN),
               rowSums(seqtab))

colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled")
rownames(track) <- sample.names

head(track)

track <- track %>% 
         data.frame() %>% 
         rownames_to_column( var = 'sample') %>%
         mutate(diff1 = filtered/input %>% round(2),
                diff2 = denoised/filtered %>% round(2),
                diff3 = merged/denoised %>% round(2),
                diff4 = tabled/merged %>% round(2),
                diff.total = tabled/input %>% round(2))

print("The median of reads kept is the following:")

summary(track$diff.total)

write_tsv(data.frame(track),str_c(output,name.run,"track_analysis.tsv"))

# The scripts generates an artifactural pdf in the main directory. 
# We are going to remove it. It is a bad solution, but a solution after all!
if (file.exists('Rplots.pdf')) file.remove('Rplots.pdf')

# END
