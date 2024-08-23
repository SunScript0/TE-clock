library(tidyverse)
library(data.table)
library(ChIPseeker)
library(GenomicFeatures)

setwd("/scratch/fmorandi/internal/RE_clock")

paths = list()
paths$data = "/scratch/fmorandi/extras/ChromAcc-clock/data/pipeline_outputs/atac"
paths$processed = "./data/processed"
paths$results = "./results/final"
paths$meta = "/scratch/fmorandi/extras/ChromAcc-clock/data/paper_data/meta_final.tsv"
# paths$gtf = "/scratch/fmorandi/external/references/GRCm39-mm39-Ensembl/Mus_musculus.GRCm39.108.gtf"

#### PLOTTING SETTINGS ####

w = 174 # mm
h = 230
w_in = w*0.0393701
h_in = h*0.0393701

#### FUNCTIONS ####

split_rep_id = function(rep_id) {
  # Split rep_id into superf, fam, id properly
  tmp = rep_id
  rep_id = str_extract(tmp, "[[:digit:]]+$")
  tmp = str_replace(tmp, "/[[:digit:]]*$", "")
  fam = str_extract(tmp, "[^/]+$")
  tmp = str_replace(tmp, "/[^/]+$", "")
  superf = tmp
  class = str_replace(tmp, "/[^/]+$", "")
  out = data.frame(
    class = class,
    superf = superf,
    fam = fam,
    rep_id = rep_id
  )
  return(out)
}

tpm = function(counts, lengths, features_are_rows=T) {
  if (features_are_rows) {
    tmp = 1e3 * sweep(counts, 1, lengths, "/")
    tmp = 1e6 * sweep(tmp, 2, colSums(tmp), "/")
  } else {
    tmp = 1e3 * sweep(counts, 2, lengths, "/")
    tmp = 1e6 * sweep(tmp, 1, rowSums(tmp), "/")
  }
  return(tmp)
}

#### LOAD DATA ####

# Read counts
# counts_ocr = fread(paste0(paths$data, "/07_peaksets_and_tables/counts_combined.tsv"),
#                    header=T, data.table=F)
counts_re = fread(paste0(paths$data, "/07_peaksets_and_tables/counts_rtes.tsv"),
                  header=T, data.table=F)
# Separate feature info
# pinfo = counts_ocr[, 1:6]
reinfo = counts_re[, 1:6]
colnames(pinfo) = c("coord", "chr", "start", "end", "strand", "length")
# counts_ocr = counts_ocr[, -c(1:6)]
counts_re = counts_re[, -c(1:6)]
# rownames(counts_ocr) = rownames(pinfo) = paste0("p", 1:nrow(pinfo))
# Tidy sample names
# all.equal(colnames(counts_ocr), colnames(counts_re))
# colnames(counts_ocr) = str_extract(colnames(counts_ocr), "/([^/]+)_clean.bam", group=1)
colnames(counts_re) = str_extract(colnames(counts_re), "/([^/]+)_clean.bam", group=1)

# Read meta
meta = read.table(paths$meta, header = T)

# Read featureCounts report
fcreport = read.table(paste0(paths$data, "/07_peaksets_and_tables/counts_rtes.tsv.summary"), header=T, row.names=1)
colnames(fcreport) = str_extract(colnames(fcreport), "\\.([^\\.]+)_clean.bam", group=1)
fcreport = data.frame(t(fcreport))
fcreport = fcreport[, colSums(fcreport)!=0]

## Read mapping logs
# logs = dir(paste0(paths$data, "/A_mapping_logs"), full.names=T)
# qc = data.frame()
# for (log in logs) {
#   f = read_file(log)
#   sname = str_extract(log, "\\/[0-9]{4}-[0-9]{2}-[0-9]{2}_([^\\/]+)\\.txt", group=1)
#   # TrimGalore! prints the total number of pairs processed after validation
#   #   This is the number of raw pairs
#   #   Sequences that were trimmed too much are removed
#   s = str_extract(f, "Total number of sequences analysed: ([0-9]+)", group=1)
#   qc[sname, "raw"] = as.numeric(s)
#   # Bowtie2 summary includes the number of pairs processed
#   #   This corresponds to the number of pairs after trimming
#   s = str_extract(f, "([0-9]+) reads; of these:", group=1)
#   qc[sname, "trimmed"] = as.numeric(s)
#   # Flagstat prints the number of properly paired reads
#   #   I can confirm that this is the number of reads seen from now on
#   #   By comparing to the remove_duplicates log
#   s = str_extract(f, "([0-9]+) \\+ [0-9]+ properly paired", group=1)
#   qc[sname, "proper_pair"] = as.numeric(s) / 2
#   # I had the pipeline print the number of mito reads from idxstats
#   s = str_extract(f, "Found ([0-9]+) mitochondrial pairs", group=1)
#   qc[sname, "mitochondrial"] = as.numeric(s)
#   # Picard MarkDuplicates prints the number of dupes in the log (divide by 2 to get pairs)
#   s = str_extract(f, "Marking ([0-9]+) records as duplicates", group=1)
#   qc[sname, "duplicates"] = as.numeric(s) / 2
#   # I had the pipeline print the number of clean reads
#   s = str_extract(f, "Clean BAM contains ([0-9]+) pairs", group=1)
#   qc[sname, "clean"] = as.numeric(s)
# }
# # Verify that qc and fc_report match up
# all(rownames(fcreport) == rownames(qc))
# all(rowSums(fcreport)/2 == qc$clean)
# 
# # Get frip
# qc$cut_sites = fcreport$Assigned + fcreport$Unassigned_NoFeatures
# qc$cut_sites_in_peak = fcreport$Assigned
# qc$frip = qc$cut_sites_in_peak / qc$cut_sites
# 
# # Combine
# meta = merge(meta, qc, by.x="SampleID", by.y=0)
# rm(qc, fcreport)
# write.table(meta, paste0(paths$results, "/qc_summary_combined.tsv"), sep="\t", quote=F, row.names=F)

#### COLLAPSE REPEATS ####

reinfo[c("class", "superf", "fam", "rep_id")] = split_rep_id(reinfo$Geneid)
counts_re[c("class", "superf", "fam", "length")] = reinfo[c("class", "superf", "fam", "Length")]
counts_re = counts_re %>%
  mutate(count = 1) %>%
  group_by(class, superf, fam) %>%
  summarise_all(sum)
reinfo = data.frame(counts_re[c("class", "superf", "fam", "count", "length")])
counts_re = counts_re %>%
  ungroup() %>%
  dplyr::select(-c(class, superf, fam, count, length))
counts_re = data.frame(counts_re)
rownames(counts_re) = rownames(reinfo) = paste0("re", 1:nrow(reinfo))
#### CHIP SEEKER ####
# 
# # ChipSeeker
# pinfo$start = as.integer(pinfo$start)
# pinfo$end = as.integer(pinfo$end)
# ranges =  GRanges(
#   seqnames = pinfo$chr,
#   ranges = IRanges(pinfo$start, pinfo$end)
# )
# txdb = makeTxDbFromGFF(paths$gtf)
# anno = annotatePeak(ranges, TxDb = txdb, tssRegion = c(-1000, 1000))
# annoDF = as.data.frame(anno) %>%
#   mutate(annotation = str_replace_all(annotation, " \\(.*\\)", "")) %>%
#   mutate(annotation = str_replace_all(annotation, "Distal ", "")) %>%
#   mutate(annotation2 = fct_recode(annotation,
#                                   "Other Genic" = "Downstream",
#                                   "Other Genic" = "3' UTR",
#                                   "Other Genic" = "5' UTR"), .after=annotation) %>%
#   dplyr::select(-(width:strand), -(geneChr:geneStrand), -transcriptId)
# pinfo = pinfo %>%
#   rownames_to_column("id") %>%
#   merge(., annoDF, by.x=c("chr", "start", "end"), by.y=c("seqnames", "start", "end"), all=T) %>%
#   column_to_rownames("id")
# pinfo = pinfo[rownames(counts_ocr), ]
# 
# all(rownames(pinfo) == rownames(counts_ocr))
# all(rownames(reinfo) == rownames(counts_re))

#### REMOVE OUTLIERS ####

all.equal(colnames(counts_re), meta$SRR_atac)
meta = subset(meta, !is.na(SRR_atac))
rownames(meta) = meta$SRR_atac
meta = meta[colnames(counts_re), ]
all.equal(colnames(counts_re), meta$SRR_atac)

meta = subset(meta, PassesQC_atac)
counts_re = counts_re[, meta$SRR_atac]


#### NORMALIZE ####

# norm_ocr = tpm(counts_ocr, pinfo$length, features_are_rows = T)
norm_re = tpm(counts_re, reinfo$length, features_are_rows = T)
reinfo$meanAcc = rowMeans(norm_re)
reinfo$meanLen = reinfo$length / reinfo$count

#### SAVE ####

atac_re_counts = counts_re
atac_re_tpm = norm_re
rinfo_atac = reinfo
meta_atac = meta
atac_re_tpm = atac_re_tmp
save(atac_re_counts, atac_re_tpm, rinfo_atac, meta_atac,
     file=paste0(paths$processed, "/atac_data.Rdata"))


