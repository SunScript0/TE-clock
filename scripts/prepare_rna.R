library(tidyverse)
library(data.table)
library(edgeR)
library(ChIPseeker)
library(ggpubr)
library(patchwork)

setwd("/scratch/fmorandi/internal/RE_clock")

paths = list()
paths$data = "./data/rna_seq"
paths$processed = "./data/processed"
paths$results = "./results/final"

#### LOAD DATA ####

# Read counts
files = dir(paste0(paths$data, "/05_counts"), pattern="*cntTable", full.names = T)
counts = list()
for (f in files) {
  sname = str_extract(f, "/([^/]+).cntTable", group=1)
  counts[[sname]] = fread(f, col.names = c("Feature", sname))
}
counts = Reduce(function(x,y) merge(x, y, by="Feature"), counts)
counts = column_to_rownames(counts, "Feature")

# Take out gene info
ginfo = data.frame(
  "Geneid" = rownames(counts), 
  "is_re" = grepl(":", rownames(counts)),
  row.names = rownames(counts))

# Keep genes exressed in at least 3 samples
counts = counts[rowSums(counts > 0) > 3, ]
ginfo = ginfo[ginfo$Geneid %in% rownames(counts), ]

# Read meta
meta = read.table(paste0(paths$data, "/meta.tsv"), sep="\t", header=T)

# Read qc info from aligment
qc = list()
files = dir(paste0(paths$data, "/04_mapped"), pattern="*Log.final.out", full.names = T)
for (f in files) {
  sname = str_extract(f, "/([^/]+)_Log.final.out", group=1)
  tmp = read_file(f)
  input_reads = str_extract(tmp, "Number of input reads \\|\\t([[:digit:]]+)", group=1)
  uniquely_mapped = str_extract(tmp, "Uniquely mapped reads number \\|\\t([[:digit:]]+)", group=1)
  multi_mapped = str_extract(tmp, "Number of reads mapped to multiple loci \\|\\t([[:digit:]]+)", group=1)
  unmapped = str_extract_all(tmp, "Number of reads unmapped.*([[:digit:]]+)", simplify = T)
  unmapped = str_extract_all(unmapped, "[[:digit:]]+", simplify = T)
  qc[[sname]] = list(
    input_reads = as.numeric(input_reads),
    uniquely_mapped = as.numeric(uniquely_mapped),
    multi_mapped = as.numeric(multi_mapped),
    unmapped = sum(as.numeric(unmapped))
  )
}

# Read qc info from counting
files = dir(paste0(paths$data, "/B_TEcounts_logs"), pattern="*txt", full.names = T)
for (f in files) {
  sname = str_extract(f, "/\\d+-\\d+-\\d+_([^/]+).txt", group=1)
  tmp = read_file(f)
  annotated = str_extract(tmp, "Total annotated reads = (\\d+) ", group=1)
  unannotated = str_extract(tmp, "Total unannotated reads = (\\d+) ", group=1)
  qc[[sname]][["annotated"]] = as.numeric(annotated)
  qc[[sname]][["unannotated"]] = as.numeric(unannotated)
}

# Merge qc and meta
qc = as.data.frame(do.call(rbind, qc)) %>%
  mutate_all(as.numeric)
meta = merge(meta, qc, by.x="SRR_rna", by.y=0)
counts = counts[, meta$SRR_rna]
write.table(meta, paste0(paths$results, "/qc_summary_rna.tsv"), sep="\t", quote=F, row.names=F)
rm(unmapped, annotated, input_reads, multi_mapped, unannotated, uniquely_mapped, qc)

#### SPLIT ginfo and re info ####

all.equal(rownames(counts), ginfo$Geneid)
counts_re = counts[ginfo$is_re, ]
counts_ge = counts[!ginfo$is_re, ]

rinfo = ginfo[ginfo$is_re, ]
ginfo = ginfo[!ginfo$is_re, ]

rm(counts)

#### CONVERT GENE IDS ####

ginfo$symbol = mapIds(org.Hs.eg.db, keys = ginfo$Geneid, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
ginfo$entrez = mapIds(org.Hs.eg.db, keys = ginfo$Geneid, keytype = "ENSEMBL", column = "ENTREZID", multiVals = "first")

# Discard genes with no symbol
ginfo = drop_na(ginfo)
ginfo = ginfo[!ginfo$symbol == "", ]
ginfo = ginfo[!duplicated(ginfo$symbol), ]
counts_ge = counts_ge[ginfo$Geneid, ]
rownames(ginfo) = ginfo$Geneid
rownames(counts_ge) = ginfo[rownames(counts_ge), "symbol"]
rownames(ginfo) = ginfo$symbol

#### SPLIT REPEAT NAMES ####

tmp = str_split_fixed(rinfo$Geneid, ":", n = 3)
rinfo$class = tmp[, 3]
rinfo$superf = paste(tmp[, 3], tmp[, 2], sep="/")
rinfo$fam = tmp[, 1]

#### NORMALIZE ####

lib_sizes = colSums(counts_ge)
cpm_ge = 1e6 * (counts_ge / lib_sizes)
cpm_re = 1e6 * (counts_re / lib_sizes)
ginfo$meanExpr = rowMeans(cpm_ge)
rinfo$meanExpr = rowMeans(cpm_re)

#### OUTLIER REMOVAL ####

meta_rna = meta

all.equal(rownames(pca$x), meta_rna$SRR_rna)
pca = prcomp(t(cpm_ge), scale=T)
pca = cbind(meta_rna, pca$x[, c("PC1", "PC2")])

ggplot(pca, aes(PC1, PC2, label = SRR_rna))+
  geom_text()
meta_rna$Outlier = meta_rna$SRR_rna %in% c(
  "SRR17465624", "SRR19068723", "SRR17465623", "SRR17465608"
)
meta_rna = subset(meta_rna, !Outlier)
pca = prcomp(t(cpm_ge[, meta_rna2$SRR_rna]), scale=T)
pca = cbind(meta_rna2, pca$x[, c("PC1", "PC2")])
ggplot(pca, aes(PC1, PC2, label = SRR_rna))+
  geom_text()

counts_ge = counts_ge[, meta_rna$SRR_rna]
counts_re = counts_re[, meta_rna$SRR_rna]
cpm_ge = cpm_ge[, meta_rna$SRR_rna]
cpm_re = cpm_re[, meta_rna$SRR_rna]

#### SAVE ####

save(counts_ge, counts_re, cpm_ge, cpm_re, ginfo, rinfo, meta_rna,
     file=paste0(paths$processed, "/rna_data.Rdata"))
