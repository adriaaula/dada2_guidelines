# DADA2 pre-processing. Q score distribution

## Load packages

library(stringr)
library(tidyverse)
library(dada2)

cat(paste0('\n',"You are using DADA2 version ", packageVersion('dada2'),'\n\n'))

cat('################################\n\n')

## Load Rscript args

args <- commandArgs(trailingOnly = TRUE)

## Variables to change for EACH RUN

path<- args[1] # where the fastq files are located
output <- args[2] # where we will dump the results//diagnostic
run.name <- args[3] # Run name. A directory for the values will be created

## Taking filenames

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

### Since a qprofile dir can be not present we check and create it

dir.create(file.path(output, "00_qprofiles"), showWarnings = FALSE)
dir.create(file.path(output, "00_qprofiles",run.name), showWarnings = FALSE)

out.diag <- file.path(output, "00_qprofiles", run.name)

## Qscore distribution

### NOT RELEVANT. Helper functions to improve  plots

lines_seq <- function(x) {
  y <-  x%%10
  if (y > 5) {
    upper <- x+(10-y)
  } else {
    upper <- x-y
  }
  lower <- upper - 50
  return(seq(upper, lower, -10))
}

improve_headers <- theme(strip.background = element_blank(),
      strip.text.x =element_text(margin = margin(.05, 0, .1, 0, "cm")),
      strip.text = element_text(size = 10),
      axis.text.y = element_text(size = 10))

### Forward

qplot_f <- plotQualityProfile(fnFs[1:9])

qplot_f_lines <- qplot_f +
                  geom_vline(xintercept = lines_seq(max(qplot_f$data$Cycle)), alpha = 0.3, size = 0.3, linetype = 'dashed') +
                  improve_headers +
                  ggtitle(paste0(run.name," - ","Forward reads"))

ggsave(plot = qplot_f_lines,
       path= out.diag,
       device="pdf",
       filename = "forward.pdf",
       width = 410,
       height = 300,
       units = 'mm')

### Same for reverse reads!

qplot_r <- plotQualityProfile(fnRs[1:9])

qplot_r_lines <- qplot_r +
                  geom_vline(xintercept = lines_seq(max(qplot_r$data$Cycle)), alpha = 0.3, size = 0.3, linetype = 'dashed') +
                  improve_headers +
                  ggtitle(paste0(run.name," - ","Reverse reads"))

ggsave(plot = qplot_r_lines,
       path= out.diag,
       device="pdf",
       filename = "reverse.pdf",
       width = 410,
       height = 300,
       units = 'mm')

cat(paste0('# Two files ("forward.pdf" and "reverse.pdf") were created in "', out.diag,'". They contain the quality profile of the first 9 samples.\n'))
cat('# You have to copy them to your local computer using "scp user@marbits:path/to/files/*pdf ." and check them to decide where to trim.\n')

cat('\nAll done!\n\n')

# The script generates an artifactual pdf in the main directory.
# We are going to remove it. It is a bad solution, but a solution after all!
if (file.exists('Rplots.pdf')) file.remove('Rplots.pdf')

# END
