library(data.table)
library(tidyverse)

setwd("/scratch/fmorandi/internal/RE_clock/")

class_order = c("LINE", "SINE", "LTR", "DNA", "DNA?", "Retroposon", "RC", 
                "tRNA", "srpRNA", "snRNA", "scRNA", "rRNA", 
                "Satellite", "Simple_repeat", "Low_complexity",
                "Unspecified", "Unknown")

define_rep_ids = function(saf_path, class_order) {
  anno = fread(saf_path, data.table = F, skip=1)
  re_names = str_replace(anno[,1], "/[^/]+$", "")
  re_names = unique(re_names)
  class = str_replace(re_names, "/.*", "")
  superf = str_replace(re_names, "/[^/]+$", "")
  fam = str_extract(re_names, "/([^/]+)$", group=1)
  re_ids = data.frame(
    class = class,
    superf = superf,
    fam = fam,
    name = re_names)
  if (!all(class  %in% class_order)) {
    stop("class order incomplete")
  }
  re_ids = re_ids %>%
    mutate(class = factor(class, levels=class_order)) %>%
    arrange(class, superf, fam)
  rownames(re_ids) = paste0("re", 1:nrow(re_ids))
  return(re_ids)
}

#### GRCh37 ####

re_ids = define_rep_ids("/scratch/fmorandi/external/references/GRCh37-hg19-Ensembl/RepeatMaskerOut/GRCh37_repeats.saf",
                        class_order)
write.csv(re_ids, "./data/annotation/GRCh37_repIds.csv")

#### GRCm38 ####

re_ids = define_rep_ids("/scratch/fmorandi/external/references/GRCm38-mm10-Ensembl/RepeatMaskerOut/GRCm38_repeats.saf",
                        class_order)
write.csv(re_ids, "./data/annotation/GRCm38_repIds.csv")
