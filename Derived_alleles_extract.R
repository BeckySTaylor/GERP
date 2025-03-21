# Control file format - a .txt file where each line is the following:
# VCF file path
# Rates file path
# Minimum GERP score
# Output file name (without file extension)
# Number of threads/cores

# Assign command-line arguments to a variable
args = commandArgs(trailingOnly=TRUE)

# Check if no arguments have been provided. If not, print message
#   that control file needs to be provided
if (length(args)==0) {
  stop("No control file provided", call.=FALSE)
}

library(vcfR)
library(tidyverse)
library(data.table)
library(magrittr)
library(parallel)

# Read in the name of the control file (given as argument 1)
# Each line of the 'ctrl' data frame is an element from the
#  control file.
ctrl = read.delim(args[1],header=F)
ctrl[,1] = as.character(ctrl[,1])

vcf <- read.vcfR(ctrl[1,1])
vcf_field_names(vcf, tag = "FORMAT")
Z <- vcfR2tidy(vcf, format_fields = c("GT", "DP"))
gt2 = Z$gt

gt2 %<>%
  filter(!str_detect(gt_GT_alleles, "\\*")) %>%
  mutate("ref" = str_remove(gt_GT_alleles, "(/|\\|).$"),
         "tar" = str_remove(gt_GT_alleles, "^.(/|\\|)"));

sca = fread(ctrl[2,1], nThread = as.numeric(ctrl[5,1])) %>%
  mutate(Pos = Pos + 1) %>%
  filter(RS_score != -1) %>%
  filter(RS_score > as.numeric(ctrl[3,1]))

sams = gt2$Indiv %>% unique();

sp1 = names(sca)[4]
sp2 = names(sca)[5]
sp3 = names(sca)[6]

out = parallel::mclapply(seq_len(length(sams)),
                         FUN = function(i){
                           system(sprintf('echo "processing %s\n"', sams[i]));

                           gt2 %>%
                             filter(Indiv == sams[i]) %>%
                             select(POS, ref, tar) -> tmp;

                           sca %>%
                             filter(Pos %in% tmp$POS) -> tmp2;

                           left_join(tmp2, tmp, by = c("Pos" = "POS")) %>%
                             filter((tar !=  !!as.symbol(sp1) &
                                       tar != !!as.symbol(sp2) &
                                       tar != !!as.symbol(sp3)) |
                                      (ref != !!as.symbol(sp1) &
                                         ref != !!as.symbol(sp2) &
                                         ref != !!as.symbol(sp3))) %>%
                             mutate("sam" = sams[i]) -> out;

                           fwrite(out, paste0(i, "_", ctrl[4,1],".txt"));#write to drive

                           return(out);
                         }, mc.preschedule = F, mc.cores = as.numeric(ctrl[5,1]))
