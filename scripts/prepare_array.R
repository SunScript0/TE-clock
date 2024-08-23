library(plyr)
library(tidyverse)
library(gtools)
library(GenomicRanges)
library(EnsDb.Hsapiens.v75)
library(ChIPseeker)
library(data.table)
library(Repitools)
library(BSgenome.Hsapiens.UCSC.hg19)
library(stringi)

setwd("/scratch/fmorandi/internal/RE_clock/")

paths = list()
paths$raw = "/scratch/fmorandi/internal/RE_clock/data/infinium"
paths$processed = "/scratch/fmorandi/internal/RE_clock/data/processed"
paths$probe_info = "/scratch/fmorandi/internal/RE_clock/data/annotation/short_probe_info_450.csv"
paths$re_info = "/scratch/fmorandi/external/references/GRCh37-hg19-Ensembl/RepeatMaskerOut/GRCh37_repeats.saf"
paths$re_ids = "/scratch/fmorandi/internal/RE_clock/data/annotation/GRCh37_repIds.csv"

##### GSE40279 #####

# Get first 100 lines
header = readLines(paste0(paths$raw, "/GSE40279/GSE40279_series_matrix.txt"), n=100)

ids = header[grepl("!Sample_geo_accession", header)]
ids = str_extract_all(ids, "GSM[0-9]+", simplify = T)

meta = header[grepl("!Sample_characteristics", header)]
meta = gsub('"', "", meta)
meta = str_split(meta, pattern = "\t", simplify = T)
meta = data.frame(t(meta))
meta = meta[-c(1), ]
colnames(meta) = c("age", "source", "plate", "sex", "ethnicity", "tissue")
meta = meta %>%
  mutate(age = str_replace_all(age, "age \\(y\\): ", "")) %>%
  mutate(source = str_replace_all(source, "source: ", "")) %>%
  mutate(plate = str_replace_all(plate, "plate: ", "")) %>%
  mutate(sex = str_replace_all(sex, "gender: ", "")) %>%
  mutate(ethnicity = str_replace_all(ethnicity, "ethnicity: ", "")) %>%
  mutate(tissue = str_replace_all(tissue, "tissue: ", "")) %>%
  mutate(age = as.numeric(age))
rownames(meta) = ids

mat_start = which(grepl("series_matrix_table_begin", header))
data = fread(paste0(paths$raw, "/GSE40279/GSE40279_series_matrix.txt"), skip = mat_start)

save(data, meta, file = paste0(paths$raw, "/GSE40279/GSE40279.Rdata"))

rm(data, ids, meta, header)

##### GSE64495 #####

# Get first 100 lines
header = readLines(paste0(paths$raw, "/GSE64495/GSE64495_series_matrix.txt"), n=100)

ids = header[grepl("!Sample_geo_accession", header)]
ids = str_extract_all(ids, "GSM[0-9]+", simplify = T)

meta = header[grepl("!Sample_characteristics", header)]
meta = gsub('"', "", meta)
meta = str_split(meta, pattern = "\t", simplify = T)
meta = data.frame(t(meta))
meta = meta[-c(1), -c(4:7, 9:10)]
colnames(meta) = c("tissue", "sex", "age", "health")
meta = meta %>%
  mutate(age = str_replace_all(age, "age: ", "")) %>%
  mutate(sex = str_replace_all(sex, "Sex: ", "")) %>%
  mutate(sex = ifelse(sex == "female", "F", "M")) %>%
  mutate(tissue = str_replace_all(tissue, "tissue: ", "")) %>%
  mutate(health = str_replace_all(health, "disease status: ", "")) %>%
  mutate(age = as.numeric(age))
rownames(meta) = ids

mat_start = which(grepl("series_matrix_table_begin", header))
data = fread(paste0(paths$raw, "/GSE64495/GSE64495_series_matrix.txt"), skip = mat_start)

save(data, meta, file = paste0(paths$raw, "/GSE64495/GSE64495.Rdata"))

rm(data, ids, meta, header)

##### GSE87648 #####

# Get all header, it doesn't contain methylation values
header = readLines(paste0(paths$raw, "/GSE87648/GSE87648_series_matrix.txt"))

ids = header[grepl("!Sample_geo_accession", header)]
ids = str_extract_all(ids, "GSM[0-9]+", simplify = T)

meta = header[grepl("!Sample_characteristics", header)] # Rows 78, 223, 337 have a problem
meta = gsub('"', "", meta)
meta = str_split(meta, pattern = "\t", simplify = T)
meta = data.frame(t(meta))
meta = meta[-c(1), -c(7)]
colnames(meta) = c("tissue", "sex", "smoking", "health", "age", "patient")
meta[c("78", "223", "337"), 4:6] = meta[c("78", "223", "337"), 3:5]
meta[c("78", "223", "337"), "smoking"] = NA
meta = meta %>%
  mutate(age = str_replace_all(age, "age: ", "")) %>%
  mutate(sex = str_replace_all(sex, "Sex: ", "")) %>%
  mutate(tissue = str_replace_all(tissue, "cell type: ", "")) %>%
  mutate(health = str_replace_all(health, "simplified_diagnosis: ", "")) %>%
  mutate(smoking = str_replace_all(smoking, "smoking status: ", "")) %>%
  mutate(patient = str_replace_all(patient, "patient_number: ", "")) %>%
  mutate(age = as.numeric(age))
rownames(meta) = ids

data = fread(paste0(paths$raw, "/GSE87648/GSE87648_betas.txt"))
old_colnames = colnames(data)[-c(1)]
data = data.frame(data)
rownames(data) = data$V1
data = data[, -c(1)]
colnames(data) = old_colnames

# Convert ids
old_ids = header[grepl("0120H", header)] # Grep one of the old ids
old_ids = str_replace_all(old_ids, " \\[Whole blood\\]", "")
old_ids = gsub('"', "", old_ids)
old_ids = str_split(old_ids, "\t", simplify=T)
ids = as.character(ids)
names(ids) = old_ids[-c(1)]
colnames(data) = ids[colnames(data)]

meta = meta[rownames(meta) %in% colnames(data), ]

save(data, meta, file = paste0(paths$raw, "/GSE87648/GSE87648.Rdata"))

rm(data, ids, old_ids, meta, header)

##### GSE147221 #####

# Get all header, it doesn't contain methylation values
header = readLines(paste0(paths$raw, "/GSE147221/GSE147221_series_matrix.txt"))

ids = header[grepl("!Sample_geo_accession", header)]
ids = str_extract_all(ids, "GSM[0-9]+", simplify = T)

meta = header[grepl("!Sample_characteristics", header)] 
# Row 5, 106, 203, 296, 389, 493, 583 have a problem
# They are technical controls, its fine to remove at the end
meta = gsub('"', "", meta)
meta = str_split(meta, pattern = "\t", simplify = T)
meta = data.frame(t(meta))
meta = meta[-c(1), -c(3)]
colnames(meta) = c("sex", "health", "age")
meta = meta %>%
  mutate(age = str_replace_all(age, "age: ", "")) %>%
  mutate(sex = str_replace_all(sex, "Sex: ", "")) %>%
  mutate(health = str_replace_all(health, "consensus_diagnosis: ", "")) %>%
  mutate(age = as.numeric(age))
rownames(meta) = ids
meta = meta[-c(4, 105, 202, 295, 388, 492, 582), ] # Age is NA

data = fread(paste0(paths$raw, "/GSE147221/GSE147221_betas.csv"))
old_colnames = colnames(data)[-c(1)]
data = data.frame(data)
rownames(data) = data$V1
data = data[, -c(1)]
colnames(data) = old_colnames
data = data %>%
  dplyr::select(-contains("Pval"))

# Convert ids
old_ids = header[grepl("200049730046_R01C01", header)] # Grep one of the old ids
old_ids = str_replace_all(old_ids, ": genomic DNA from whole Blood", "")
old_ids = gsub('"', "", old_ids)
old_ids = str_split(old_ids, "\t", simplify=T)
ids = as.character(ids)
names(ids) = old_ids[-c(1)]
colnames(data) = ids[colnames(data)]

meta = meta[rownames(meta) %in% colnames(data), ]

save(data, meta, file = paste0(paths$raw, "/GSE147221/GSE147221.Rdata"))

rm(data, ids, old_ids, meta, header)

##### GSE157131 #####

# === Infinium meta ===
# Get all header, it doesn't contain methylation values
header = readLines(paste0(paths$raw, "/GSE157131/GSE157131-GPL13534_series_matrix.txt"))

ids = header[grepl("!Sample_geo_accession", header)]
ids = as.character(str_extract_all(ids, "GSM[0-9]+", simplify = T))
ids2 = header[grepl("!Sample_description", header)]
ids2 = as.character(str_extract_all(ids2, "AAGENOA_LCL_[0-9]+", simplify = T))

meta1 = header[grepl("!Sample_characteristics", header)]
meta1 = gsub('"', "", meta1)
meta1 = str_split(meta1, pattern = "\t", simplify = T)
meta1 = data.frame(t(meta1))
meta1 = meta1[-c(1), ]
colnames(meta1) = c("tissue", "sex", "age", "sysp", "diasp")
meta1 = meta1 %>%
  mutate(age = str_replace_all(age, "age\\(yrs\\):", "")) %>%
  mutate(sex = str_replace_all(sex, "gender: ", "")) %>%
  mutate(sysp = str_replace_all(sysp, "systolic blood pressure \\(mmhg\\): ", "")) %>%
  mutate(diasp = str_replace_all(diasp, "diastolic blood pressure \\(mmhg\\): ", "")) %>%
  mutate(tissue = str_replace_all(tissue, "tissue: ", "")) %>%
  mutate_at(c("age", "sysp", "diasp"), as.numeric)
rownames(meta1) = ids
meta1$alt_id = ids2

# === Infinium meta ===
# Get all header, it doesn't contain methylation values
header = readLines(paste0(paths$raw, "/GSE157131/GSE157131-GPL21145_series_matrix.txt"))

ids = header[grepl("!Sample_geo_accession", header)]
ids = as.character(str_extract_all(ids, "GSM[0-9]+", simplify = T))
ids2 = header[grepl("!Sample_description", header)]
ids2 = as.character(str_extract_all(ids2, "AAGENOA_LCL_[0-9]+", simplify = T))

meta2 = header[grepl("!Sample_characteristics", header)]
meta2 = gsub('"', "", meta2)
meta2 = str_split(meta2, pattern = "\t", simplify = T)
meta2 = data.frame(t(meta2))
meta2 = meta2[-c(1), ]
colnames(meta2) = c("tissue", "sex", "age", "sysp", "diasp")
meta2 = meta2 %>%
  mutate(age = str_replace_all(age, "age\\(yrs\\):", "")) %>%
  mutate(sex = str_replace_all(sex, "gender: ", "")) %>%
  mutate(sysp = str_replace_all(sysp, "systolic blood pressure \\(mmhg\\): ", "")) %>%
  mutate(diasp = str_replace_all(diasp, "diastolic blood pressure \\(mmhg\\): ", "")) %>%
  mutate(tissue = str_replace_all(tissue, "tissue: ", "")) %>%
  mutate_at(c("age", "sysp", "diasp"), as.numeric)
rownames(meta2) = ids
meta2$alt_id = ids2

# === Concat meta ===
meta = rbind(meta1, meta2)
rm(meta1, meta2)

# Define health column based on pressure
meta$health = "Hypertension stage 3"
meta[meta$sysp < 180 & meta$diasp < 120, "health"] = "Hypertension stage 2"
meta[meta$sysp < 140 & meta$diasp < 90, "health"] = "Hypertension stage 1"
meta[meta$sysp < 140 & meta$diasp < 90, "health"] = "Hypertension stage 1"
meta[meta$sysp < 130 & meta$diasp < 80, "health"] = "Elevated pressure"
meta[meta$sysp < 120 & meta$diasp < 80, "health"] = "Healthy"
meta[meta$sysp < 105 & meta$diasp < 60, "health"] = "Low"

data = fread(paste0(paths$raw, "/GSE157131/GSE157131_Matrix_processed_beta_geo_08252020.txt"), data.table=F)
data = data[, !grepl("Pval", colnames(data))]
id_conv = c("ID_REF", rownames(meta))
names(id_conv) = c("ID_REF", meta$alt_id)
colnames(data) = id_conv[colnames(data)]

save(data, meta, file = paste0(paths$raw, "/GSE157131/GSE157131.Rdata"))

rm(data, id_conv, meta, header, ids, ids2)

##### GSE51032 #####

# Get all header, it doesn't contain methylation values
header = readLines(paste0(paths$raw, "/GSE51032/GSE51032_series_matrix.txt"), n=100)

ids = header[grepl("!Sample_geo_accession", header)]
ids = str_extract_all(ids, "GSM[0-9]+", simplify = T)

meta = header[grepl("!Sample_characteristics", header)] 
meta = gsub('"', "", meta)
meta = str_split(meta, pattern = "\t")
meta = data.frame(
  sex = meta[[1]],
  age = meta[[2]],
  age_at_menarche = meta[[3]],
  cancer_type = meta[[4]],
  time_to_diagnosis = meta[[5]]
)
meta = meta[-c(1),]

all(grepl("gender", meta$sex))
all(grepl("age", meta$age))
unique(meta[grepl("cancer", meta$age_at_menarche), "time_to_diagnosis"])
meta[grepl("cancer", meta$age_at_menarche), 4:5] = meta[grepl("cancer", meta$age_at_menarche), 3:4]
meta[grepl("cancer", meta$age_at_menarche), 3] = ""
all(grepl("menarche", meta$age_at_menarche) | meta$age_at_menarche == "")
all(grepl("cancer", meta$cancer_type) | meta$cancer_type == "")
all(grepl("time", meta$time_to_diagnosis) | meta$time_to_diagnosis == "")

meta = meta %>%
  dplyr::select(-age_at_menarche) %>%
  mutate(age = str_replace_all(age, "age: ", "")) %>%
  mutate(sex = str_replace_all(sex, "gender: ", "")) %>%
  mutate(cancer_type = str_replace_all(cancer_type, "cancer type \\(icd-10\\): ", "")) %>%
  mutate(time_to_diagnosis = str_replace_all(time_to_diagnosis, "time to diagnosis: ", "")) %>%
  mutate(age = as.numeric(age)) %>%
  mutate(time_to_diagnosis = as.numeric(time_to_diagnosis))
meta[meta$cancer_type == "", "cancer_type"] = NA
meta$healthy = is.na(meta$cancer_type)
rownames(meta) = ids

mat_start = which(grepl("series_matrix_table_begin", header))
data = fread(paste0(paths$raw, "/GSE51032/GSE51032_series_matrix.txt"), skip = mat_start)
data = data.frame(data)
rownames(data) = data$ID_REF
data = data[, -c(1)]

all(rownames(meta) %in% colnames(data))
meta = meta[colnames(data), ]

save(data, meta, file = paste0(paths$raw, "/GSE51032/GSE51032.Rdata"))

rm(data, ids, meta, header)

##### CONCATENATE #####

load(paste0(paths$raw, "/GSE40279/GSE40279.Rdata"))
data = column_to_rownames(data, var = "ID_REF")
data_cat = data
meta$study = "GSE40279"
meta$health = "Healthy"
meta$ID = rownames(meta)
meta_cat = meta[c("ID", "study", "tissue", "age", "sex", "health", "ethnicity")]

load(paste0(paths$raw, "/GSE64495/GSE64495.Rdata"))
data = column_to_rownames(data, var = "ID_REF")
common_cpgs = intersect(rownames(data_cat), rownames(data))
data_cat = cbind(data_cat[common_cpgs, ], data[common_cpgs, ])
meta$study = "GSE64495"
meta$ID = rownames(meta)
meta_cat = rbind.fill(meta_cat, meta[c("ID", "study", "tissue", "age", "sex", "health")])

load(paste0(paths$raw, "/GSE87648/GSE87648.Rdata"))
common_cpgs = intersect(rownames(data_cat), rownames(data))
data_cat = cbind(data_cat[common_cpgs, ], data[common_cpgs, ])
meta$study = "GSE87648"
meta$ID = rownames(meta)
meta_cat = rbind.fill(meta_cat, meta[c("ID", "study", "tissue", "age", "sex", "health", "smoking")])

load(paste0(paths$raw, "/GSE147221/GSE147221.Rdata"))
common_cpgs = intersect(rownames(data_cat), rownames(data))
data_cat = cbind(data_cat[common_cpgs, ], data[common_cpgs, ])
meta$study = "GSE147221"
meta$tissue = "whole blood"
meta$ID = rownames(meta)
meta_cat = rbind.fill(meta_cat, meta[c("ID", "study", "tissue", "age", "sex", "health")])

load(paste0(paths$raw, "/GSE157131/GSE157131.Rdata"))
data = column_to_rownames(data, var = "ID_REF")
common_cpgs = intersect(rownames(data_cat), rownames(data))
data_cat = cbind(data_cat[common_cpgs, ], data[common_cpgs, ])
meta$study = "GSE157131"
meta$ethnicity = "African American"
meta$ID = rownames(meta)
meta_cat = rbind.fill(meta_cat, meta[c("ID", "study", "tissue", "age", "sex", "health", "ethnicity")])

# load(paste0(paths$raw, "/GSE51032/GSE51032.Rdata"))
# common_cpgs = intersect(rownames(data_cat), rownames(data))
# data_cat = cbind(data_cat[common_cpgs, ], data[common_cpgs, ])
# meta$study = "GSE51032"
# meta$tissue = "peripheral blood leukocytes"
# meta$ethnicity = NA
# meta$ID = rownames(meta)
# meta$health = meta$cancer_type
# meta_cat = rbind.fill(meta_cat, meta[c("ID", "study", "tissue", "age", "sex", "health", "ethnicity")])

data = data_cat %>%
  mutate(across(where(is.character), ~na_if(., "NULL"))) %>%
  mutate_if(is.character, as.numeric)
data = data.table::transpose(data)
rownames(data) = colnames(data_cat)
colnames(data) = rownames(data_cat)

# na_samples = colSums(is.na(data))
# na_features = rowSums(is.na(data))
# sum(na_features > 10000)

data = data[, colSums(is.na(data)) == 0]
meta = meta_cat %>%
  mutate(health = str_replace_all(health, "Control", "Healthy")) %>%
  mutate(health = str_replace_all(health, "HL", "Healthy")) %>%
  mutate(health = str_replace_all(health, "UC", "Ulcerative Colitis")) %>%
  mutate(health = str_replace_all(health, "CD", "Crohn's Disease")) %>%
  mutate(health = str_replace_all(health, "HS", "Symptomatic (Gut)")) %>%
  mutate(tissue = str_replace_all(tissue, "whole blood", "Whole blood")) %>%
  mutate(smoking = tidyr::replace_na(smoking, "Don't know"))

data = data[meta$ID, ]
all(rownames(data) == meta$ID)

rm(data_cat, meta_cat)

# <Output>
saveRDS(data, file = paste0(paths$processed, "/infinium_data.rds"))
write.table(meta, paste0(paths$processed, "/infinium_meta.tsv"), sep="\t", quote=F)

##### GENOMIC POSITION OF ALL RES #####

# Read RE annotations
anno_re = fread(paths$re_info, col.names=c("rep_id", "chr", "start", "end", "strand"),
                skip=1, data.table=F)
anno_re$fragment_id = paste0(anno_re$rep_id, ".", 1:nrow(anno_re))

# Restrict to main chromosomes
main_chr = c(1:22, "X", "Y")
anno_re = subset(anno_re, chr %in% main_chr)

# RE ranges object
RE_ranges = GRanges(
  seqnames = anno_re$chr,
  ranges = IRanges(anno_re$start, end = anno_re$end, names = anno_re$fragment_id))

# Genic position of REs in general
txdb = EnsDb.Hsapiens.v75
anno_chipseeker = annotatePeak(RE_ranges, TxDb = txdb, tssRegion = c(-1000, 1000))

# Tidy up
anno_chipseekerDF = as.data.frame(anno_chipseeker) %>%
  mutate(annotation = str_replace_all(annotation, " \\(.*\\)", "")) %>%
  mutate(annotation = str_replace_all(annotation, "Distal ", "")) %>%
  mutate(annotation2 = fct_recode(annotation, 
                                  "Other Genic" = "Downstream", 
                                  "Other Genic" = "3' UTR", 
                                  "Other Genic" = "5' UTR"), .after=annotation)
anno_chipseekerDF = anno_chipseekerDF[, c("seqnames", "start", "end", 
                                          "annotation", "annotation2", 
                                          "transcriptBiotype", "distanceToTSS", "geneId")]
anno_chipseekerDF = merge(anno_re, anno_chipseekerDF, by.x=c("chr", "start", "end"), by.y=c("seqnames", "start", "end"))

# Count for each type of location
allRE_locs = anno_chipseekerDF %>%
  group_by(annotation) %>%
  dplyr::summarize(n = n())
ggplot(allRE_locs, aes(x="", y=n, fill=annotation))+
  geom_bar(stat="identity", width=1)+
  coord_polar("y", start=0)

# <Output>
saveRDS(anno_chipseekerDF, file = paste0(paths$processed, "/re_chipseeker_GRCh37.rds"))
write.table(allRE_locs, file = paste0(paths$processed, "/re_locations.tsv"))

##### PROBE INFO #####

# Read probe information
probe_info = read.csv(paths$probe_info)
# <!> Removing any probe not in genome build 37
probe_info = subset(probe_info, probe_info$Genome_Build == 37)
# Remove uninteresting columns for clarity
probe_info = probe_info[, c(1, 12, 13, 10, 26)] #, 22:33

# Map probes to RE elements
probes = GRanges(
  seqnames = probe_info$CHR,
  ranges = IRanges(probe_info$MAPINFO, names = probe_info$IlmnID)
)
mappings = data.frame(findOverlaps(RE_ranges, probes))
mappings$queryHits=names(RE_ranges)[mappings$queryHits]
mappings$subjectHits=names(probes)[mappings$subjectHits]

# <!> I discard multiple mappings in this process (very few)
mappings = merge(mappings, anno_re, by.x="queryHits", by.y="fragment_id")
colnames(mappings) = c("fragment_name", "probe", "rep_name", "chr", "re_start", "re_end", "re_strand")

# Split re_id into superf, fam, id properly
tmp = mappings$rep_name
mappings$rep_id = str_extract(tmp, "[[:digit:]]+$")
tmp = str_replace(tmp, "/[[:digit:]]+$", "")
mappings$fam = str_extract(tmp, "[^/]+$")
tmp = str_replace(tmp, "/[^/]+$", "")
mappings$superf = tmp
mappings$class = str_replace(mappings$superf, "/[^/]+$", "")

mappings = distinct(mappings, probe, .keep_all = T)
  
probe_info = merge(probe_info, mappings, by.x="IlmnID", by.y="probe", all=T)

# ChipSeeker
txdb = EnsDb.Hsapiens.v75
anno = annotatePeak(probes, TxDb = txdb, tssRegion = c(-1000, 1000))

# Tidy up
annoDF = as.data.frame(anno) %>%
  mutate(annotation = str_replace_all(annotation, " \\(.*\\)", "")) %>%
  mutate(annotation = str_replace_all(annotation, "Distal ", "")) %>%
  mutate(annotation2 = fct_recode(annotation, 
                                  "Other Genic" = "Downstream", 
                                  "Other Genic" = "3' UTR", 
                                  "Other Genic" = "5' UTR"), .after=annotation)%>%
  mutate(annotation3 = fct_recode(annotation2,
                                  "Genic" = "Promoter",
                                  "Genic" = "Intron",
                                  "Genic" = "Exon",
                                  "Genic" = "Other Genic"))
annoDF = annoDF[, c("seqnames", "start", "annotation", "annotation2", "annotation3", "transcriptBiotype", "distanceToTSS", "geneId")]
probe_info = merge(probe_info, annoDF, by.x=c("CHR", "MAPINFO"), by.y=c("seqnames", "start"))

# Calculate CpG density around probes
probes = GRanges(
  seqnames = paste0("chr", probe_info$CHR), # Only autosomes and X,Y so this works
  ranges = IRanges(start=probe_info$MAPINFO)
)
probe_info$CpG_density = cpgDensityCalc(probes, window = 100, organism=BSgenome.Hsapiens.UCSC.hg19) / 100

# Last quick modifications
probe_info = dplyr::mutate(probe_info, re_len = re_end-re_start, .after=re_end)

# <Output>
saveRDS(probe_info, file = paste0(paths$processed, "/infinium_probes.rds"))

##### COLLAPSE REPEATS #####

data_infi = readRDS(paste0(paths$processed, "/infinium_data.rds"))
probe_info = readRDS(paste0(paths$processed, "/infinium_probes.rds"))
fam_ids = read.csv("/scratch/fmorandi/internal/RE_clock/data/annotation/GRCh37_repIds.csv", row.names = 1) %>%
  rownames_to_column("fam_id")

probe_info = merge(probe_info, fam_ids, by=c("class", "superf", "fam"), all.x=T)
rownames(probe_info) = probe_info$IlmnID
probe_info = probe_info[colnames(data_infi), ]

collapse_to_clusters = function(data, clusts) {
  uniq_clusts = sort(setdiff(unique(clusts), NA))
  data_clust = list()
  for (c in uniq_clusts) {
    cols = which(!is.na(clusts) & clusts == c)
    data_clust[[c]] = rowMeans(data[cols], na.rm=T)
  }
  data_clust = do.call(cbind, data_clust)
  rownames(data_clust) = rownames(data)
  colnames(data_clust) = uniq_clusts
  if (is.data.frame(data)) {
    data_clust = data.frame(data_clust)
  }
  return(data_clust)
}

###### Real ######

re_info1 = probe_info %>%
  mutate(clust = fam_id)
data_re1 = collapse_to_clusters(data_infi, re_info1$clust)
re_info1 = re_info1 %>%
  group_by(class, superf, fam, fam_id, clust) %>%
  dplyr::filter(!is.na(clust)) %>%
  dplyr::summarize(n=n())%>%
  column_to_rownames("clust")

re_info2 = probe_info %>%
  mutate(clust = paste(fam_id, annotation3, sep="_"))
re_info2[is.na(re_info2$class), "clust"] = NA
data_re2 = collapse_to_clusters(data_infi, re_info2$clust)
re_info2 = re_info2 %>%
  group_by(class, superf, fam, fam_id, annotation3, clust) %>%
  dplyr::filter(!is.na(clust)) %>%
  dplyr::summarize(n=n())%>%
  column_to_rownames("clust")

re_info3 = probe_info %>%
  mutate(clust = paste(fam_id, annotation2, sep="_"))
re_info3[is.na(re_info3$class), "clust"] = NA
data_re3 = collapse_to_clusters(data_infi, re_info3$clust)
re_info3 = re_info3 %>%
  group_by(class, superf, fam, fam_id, annotation2, clust) %>%
  dplyr::filter(!is.na(clust)) %>%
  dplyr::summarize(n=n())%>%
  column_to_rownames("clust")

# <Output>
save(data_re1, re_info1, file = paste0(paths$processed, "/infinium_data_collapsed1.Rdata"))
save(data_re2, re_info2, file = paste0(paths$processed, "/infinium_data_collapsed2.Rdata"))
save(data_re3, re_info3, file = paste0(paths$processed, "/infinium_data_collapsed3.Rdata"))
