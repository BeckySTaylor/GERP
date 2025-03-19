library(vcfR)
library(tidyverse)
library(data.table)
library(magrittr)
library(parallel)

vcf <- read.vcfR("./Input_VCF_Scaffold1.vcf.gz")
vcf_field_names(vcf, tag = "FORMAT")
Z <- vcfR2tidy(vcf, format_fields = c("GT", "DP"))
gt2 = Z$gt

gt2 %<>%
  filter(!str_detect(gt_GT_alleles, "\\*")) %>%
  mutate("ref" = str_remove(gt_GT_alleles, "(/|\\|).$"),
         "tar" = str_remove(gt_GT_alleles, "^.(/|\\|)"));

sca = fread("./Scaffold15.fasta.rates", nThread = 8) %>%
  mutate(Pos = Pos + 1) %>%
  filter(V2 > 0)

sams = gt2$Indiv %>% unique();

out = parallel::mclapply(seq_len(length(sams)),
                         FUN = function(i){
                           system(sprintf('echo "processing %s\n"', sams[i]));

                           gt2 %>%
                             filter(Indiv == sams[i]) %>%
                             select(POS, ref, tar) -> tmp;

                           sca %>%
                             filter(Pos %in% tmp$POS) -> tmp2;

                           left_join(tmp2, tmp, by = c("Pos" = "POS")) %>%
                             filter((tar !=  Odocoileus_virginianus &
                                      tar != Alces_alces &
                                      tar != Cervus_elaphus) |
                                      (ref !=  Odocoileus_virginianus &
                                      ref != Alces_alces &
                                      ref != Cervus_elaphus)) %>%
                             mutate("sam" = sams[i]) -> out;

                           fwrite(out, paste0(i, "_Scaffold15_over2.txt"));#write to drive

                           return(out);
                         }, mc.preschedule = F, mc.cores = 8)
