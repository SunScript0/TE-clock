library(tidyverse)
library(data.table)

setwd("/scratch/fmorandi/internal/RE_clock/")

paths = list()
paths$data = "./data/whi"
paths$out = "./data/processed"

##### META #####

# DEATHDY, ENDFOLLOWDY
meta = list()
meta[[1]] = fread(paste0(paths$data, "/meta/dem_ctos_inv.dat"), data.table=F)
meta[[2]] = fread(paste0(paths$data, "/meta/outc_aging_ctos_inv.dat"), data.table=F)
meta[[3]] = fread(paste0(paths$data, "/meta/outc_ct_os_inv.dat"), data.table=F)
meta[[4]] = fread(paste0(paths$data, "/meta/f41_imputed_ctos_inv.dat"), data.table=F)

##### BA23 #####

data = fread(paste0(paths$data, "/BA23/BA23Methylation.txt"), data.table=F)

ids1 = unlist(data[1, ])
ids2 = unlist(data[2, ])
ids = ids1
ids[is.na(ids1)] = ids2[is.na(ids1)]

length(unique(ids1))
sum(is.na(ids1))
length(ids1) - sum(is.na(ids1))
length(unique(ids))

tmp = data.frame(ids = ids[-c(1)]) %>%
  group_by(SubjectID = ids) %>%
  mutate(Replicate = dplyr::row_number()) %>%
  mutate(SampleID = paste(SubjectID, Replicate, sep = ".")) %>%
  dplyr::select(SampleID, SubjectID)
meta_whi = Reduce(function(x, y) merge(x, y, by="ID"), meta)
meta_whi = merge(tmp, meta_whi, by.x="SubjectID", by.y="ID")
  
data = data[-c(1:2), ]
colnames(data) = c("IlmnID", tmp$SampleID)

##### SAVE #####

data_whi = data

# <Output>
save(data_whi, file = paste0(paths$out, "/whi_data.Rdata"))
write.table(meta_whi, file = paste0(paths$out, "/whi_meta.tsv"))
