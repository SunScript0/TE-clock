library(methylKit)
library(tidyverse)
library(GenomicRanges)
library(data.table)
library(unixtools)
library(Repitools)
library(ChIPseeker)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)

setwd("/scratch/fmorandi/internal/RE_clock")
set.tempdir("/scratch/fmorandi/tmp/rtmp")
source("./scripts/utils.R")

paths = list()
paths$re_info = "/scratch/fmorandi/external/references/GRCm38-mm10-Ensembl/RepeatMaskerOut/GRCm38_repeats.saf"
paths$re_ids = "/scratch/fmorandi/internal/RE_clock/data/annotation/GRCm38_repIds.csv"
paths$out = "./data/processed"
paths$bis = "./data/rrbs_mouse/thompson"
paths$clocks = "./data/published_clocks"

chr_aliases = read.table("/scratch/fmorandi/external/references/GRCm38-mm10-Ensembl/mm10.chrom.aliases",
                         sep = "\t", fill=T, col.names = c("ucsc", "grc", "idk", "ens"))
rownames(chr_aliases) = chr_aliases$grc

#### COLLAPSE INDIVIDUAL DATASETS TO RE ####

##### Load methylation data #####

# As reminder: GRCm38 = mm10

# Bismark coverage (epi_memory, mat_nutrition, stubbs)
# files = Sys.glob(paste0(paths$bis, "/05_bismark_outputs/*cov.gz")) # If I did the processing
files = Sys.glob(paste0(paths$bis, "/bismark_coverage/*cov.gz")) # If downloaded pre-processed
# files = files[1:70] # Sometimes too many files to process at once
files = files[71:length(files)]
fnames = str_extract(files, "\\/([^\\/]+$)", group = 1)
fnames = str_replace_all(fnames, ".cov.gz", "")
fnames = str_replace_all(fnames, ".gz.bismark", "")
fnames = str_replace_all(fnames, "_pe$", "")
metobj = methRead(
  as.list(files),
  sample.id=as.list(fnames),
  treatment=rep(0, length(files)),
  pipeline = "bismarkCoverage",
  assembly="mm10",
  mincov = 1,
)

# Bismark overlap (petkovich)
files = Sys.glob(paste0(paths$bis, "/overlap/*overlap.txt.gz"))
fnames = str_extract(files, "\\/([^\\/]+$)", group = 1)
fnames = str_replace_all(fnames, ".overlap.txt.gz", "")
flist = list()
for (i in 1:length(files)) {
  tmp = fread(files[i], data.table=F)
  info = str_split(tmp$index, "\\|", n=5, simplify = T)
  start = as.integer(sub(":", "", info[, 5]))
  cov = as.integer(tmp[, 3])
  numCs = as.integer(round(tmp[, 3] * tmp[, 2] * 0.01))
  dt = data.frame(
    chr = chr_aliases[info[, 4], "ens"],
    start = start,
    end = start,
    strand = rep("*", length(start)),
    coverage = cov,
    numCs = numCs,
    numTs = cov - numCs
  )
  dt = dt[dt$chr != "", ]
  obj = new("methylRaw", dt, sample.id=fnames[i], assembly="mm10",
            context="CpG",resolution="base")
  flist[[i]] = obj
}
metobj = new("methylRawList", flist, treatment=rep(0, length(files)))

# BSMAP methratio (huntington)
files = Sys.glob(paste0(paths$bis, "/methratio/*methratio.txt.gz"))
fnames = str_extract(files, "\\/([^\\/]+$)", group = 1)
fnames = str_replace_all(fnames, ".methratio.txt.gz", "")
fnames = str_extract(fnames, "^([[:alnum:]]+)_", group=1)
pipeline = list(
  fraction = T,
  chr.col = 1,
  start.col = 2,
  end.col=2,
  coverage.col = 6,
  strand.col = 3,
  freqC.col = 5
)
metobj = methRead(
  as.list(files),
  sample.id=as.list(fnames),
  treatment=rep(0, length(files)),
  pipeline = pipeline,
  assembly="mm10",
  mincov = 1
)

# BSSEEKER CGfinal (thompson)
files = Sys.glob(paste0(paths$bis, "/cgfinal/*CGfinal.txt.gz"))
files = files[grepl("muscle", files, ignore.case = T)]
fnames = str_extract(files, "\\/([^\\/]+$)", group = 1)
fnames = str_replace_all(fnames, ".CGfinal.txt.gz", "")
fnames = str_extract(fnames, "^([[:alnum:]]+)_", group=1)
flist = list()
for (i in 1:length(files)) {
  print(round(100*i/length(files)))
  tmp = fread(files[i], data.table=F)
  dt = data.frame(
    chr = tmp[, 1],
    start = tmp[, 3],
    end = tmp[, 3],
    strand = rep("*", nrow(tmp)),
    coverage = tmp[, 5],
    numCs = tmp[, 4],
    numTs = tmp[, 5] - tmp[, 4]
  )
  obj = new("methylRaw", dt, sample.id=fnames[i], assembly="mm10",
            context="CpG",resolution="base")
  flist[[i]] = obj
}
metobj = new("methylRawList", flist, treatment=rep(0, length(files)))

##### Convert to Ensembl chromosome names #####

# Needed for: huntington, thompson
tmp = unique(metobj[[1]]$chr)
chr_conv = str_replace_all(tmp, "chrM", "chrMT")
chr_conv = str_replace_all(chr_conv, "chr", "")
names(chr_conv) = tmp
for (i in 1:length(metobj)) {
  metobj[[i]]$chr = chr_conv[metobj[[i]]$chr]
  metobj[[i]] = metobj[[i]][!is.na(metobj[[i]]$chr), ]
}

##### Get individual CpGs for published clocks #####

clocks = list()
clocks$thompson = read.table(paste0(paths$clocks, "/thompson_coefs.txt"), header=T, sep="\t")
clocks$meer = read.table(paste0(paths$clocks, "/meer_coefs.txt"), header=T, sep="\t")
clocks$meer$Chr = str_replace(clocks$meer$Chr, "chr", "")
clock_cpgs = do.call(rbind, clocks)
clock_cpgs = subset(clock_cpgs, Chr != "(Intercept)")
clock_cpgs = distinct(clock_cpgs, Chr, Coord)

pubclock_ranges = GRanges(
  seqnames = clock_cpgs$Chr,
  ranges = IRanges(
    start = clock_cpgs$Coord,
    end=clock_cpgs$Coord))

# Select CpGs overlappingREs and unite
metobj_clocks = selectByOverlap(metobj, pubclock_ranges)
getCoverageStats(metobj[[1]], plot=TRUE)
getCoverageStats(metobj_clocks[[1]], plot=TRUE)
col_table = methylKit::unite(metobj_clocks, min.per.group = 1L)
all(attr(col_table, "sample.ids") == fnames)
col_table = getData(col_table) %>%
  dplyr::select(-strand, -starts_with("coverage"))

new_names = paste0(rep(c("numCs_", "numTs_"), length(fnames)), rep(fnames, each=2))
new_names = c("chr", "start", "end", new_names)
colnames(col_table) = new_names

saveRDS(col_table, paste0(paths$bis, "/single_muscle.rds"))

##### Load repeat annotations #####

anno_re = fread(paths$re_info, col.names=c("rep_id", "chr", "start", "end", "strand"),
                skip=1, data.table=F)
anno_re = distinct(anno_re, chr, start, end, .keep_all=T)
anno_re$fragment_id = paste0(anno_re$rep_id, ".", 1:nrow(anno_re))
fam_ids = read.csv(paths$re_ids, row.names = 1) %>%
  rownames_to_column("fam_id")

##### Collapsing to fragment #####

# Useful methods
# getMethylationStats(metobj[[1]], plot=TRUE)
# getCoverageStats(metobj[[1]], plot=TRUE)
# metobj2 = selectByOverlap(metobj, anno_ranges)

# getCoverageStats(metobj[[1]], plot=TRUE)

# First calculate % methylation, then average, since there may be amplification bias for different regions
# Not implemented yet <!>
avgMets = sapply(metobj, function(x) mean(x$numCs / (x$numCs + x$numTs)))

anno_ranges = GRanges(
  seqnames = anno_re$chr,
  names = anno_re$fragment_id,
  ranges = IRanges(
    start = anno_re$start,
    end = anno_re$end))

# Select CpGs overlappingREs and collapse
metobj2 = selectByOverlap(metobj, anno_ranges)
collapsed = regionCounts(metobj2, anno_ranges, destrand=T)
# getCoverageStats(collapsed[[1]], plot=TRUE)

col_table = methylKit::unite(collapsed, min.per.group = 1L)
all(attr(col_table, "sample.ids") == fnames)
col_table = getData(col_table)
col_table = anno_re %>%
  dplyr::select(-strand) %>%
  merge(., col_table, by=c("chr", "start", "end")) %>%
  dplyr::select(-starts_with("coverage"))

new_names = paste0(rep(c("numCs_", "numTs_"), length(fnames)), rep(fnames, each=2))
new_names = c("chr", "start", "end", "rep_id", "fragment_id", "strand", new_names)
colnames(col_table) = new_names

saveRDS(col_table, paste0(paths$bis, "/individual_kidney.rds"))

##### Chipseeker #####

# Restrict to main chromosomes
main_chr = c(1:22, "X", "Y")
anno_re = subset(anno_re, chr %in% main_chr)

# RE ranges object
anno_ranges = GRanges(
  seqnames = anno_re$chr,
  names = anno_re$fragment_id,
  ranges = IRanges(
    start = anno_re$start,
    end = anno_re$end))

# Genic position of REs in general
txdb = EnsDb.Mmusculus.v79
anno_chipseeker = annotatePeak(anno_ranges, TxDb = txdb, tssRegion = c(-1000, 1000))

# Tidy up
anno_chipseekerDF = as.data.frame(anno_chipseeker) %>%
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
anno_chipseekerDF = anno_chipseekerDF[, c("seqnames", "start", "end", 
                                          "annotation", "annotation2", "annotation3",
                                          "transcriptBiotype", "distanceToTSS", "geneId")]
anno_chipseekerDF = merge(anno_re, anno_chipseekerDF, by.x=c("chr", "start", "end"), by.y=c("seqnames", "start", "end"))

saveRDS(anno_chipseekerDF, file = paste0(paths$out, "/re_chipseeker_GRCm38.rds"))

##### Collapsing to family x region #####

anno_chipseekerDF[c("superf", "fam", "rep_id")] = split_rep_id(anno_chipseekerDF$rep_id)
anno_chipseekerDF = merge(anno_chipseekerDF, fam_ids, by=c("superf", "fam"))
rownames(anno_chipseekerDF) = anno_chipseekerDF$fragment_id
anno_chipseekerDF$clust = paste(anno_chipseekerDF$fam_id, anno_chipseekerDF$annotation3, sep="_")

# col_table = readRDS("./data/rrbs_mouse/thompson/individual_muscle.rds")

col_table$clust = anno_chipseekerDF[col_table$fragment_id, "clust"]
rep_table = col_table %>%
  dplyr::select(-chr, -start, -end, -rep_id, -fragment_id, -strand) %>%
  group_by(clust) %>%
  summarize_all(sum, na.rm=T)

saveRDS(rep_table, "./data/rrbs_mouse/thompson/collapsed_muscle.rds")

##### RE bed file #####

# bed = data.frame(
#   chr = anno$chr,
#   start = anno$start,
#   end = anno$end,
#   name = anno$repID
# )
# 
# gz1 = gzfile(paste0(paths$out, "GRCm38_REs.bed.gz"), "w")
# write.table(bed, gz1, quote=F, sep="\t", row.names=F)
# close(gz1)

##### Save RE genomic stats #####

anno2 = anno_re %>%
  mutate(chr_ucsc = paste0("chr", chr)) %>%
  mutate(chr_ucsc = str_replace_all(chr_ucsc, "chrMT", "chrM"))  %>%
  dplyr::filter(!grepl("\\.", chr_ucsc)) # Canonical chr only
re_ranges = GRanges(
  seqnames = anno2$chr_ucsc,
  ranges = IRanges(start=anno2$start, end=anno2$end)
)
anno2$CpGN = cpgDensityCalc(re_ranges, organism=BSgenome.Mmusculus.UCSC.mm10)

anno2[c("superf", "fam", "rep_id")] = split_rep_id(anno2$rep_id)
anno2 = anno2 %>%
  group_by(superf, fam) %>%
  summarise(n = n(), ncpg = sum(CpGN))
anno3 = merge(anno2, fam_ids, by=c("superf", "fam")) %>%
  dplyr::select(fam_id, class, superf, fam, name, n, ncpg)
  
# <Output>
write.table(anno3, paste0(paths$out, "/re_cpg_counts.tsv"), sep="\t", quote=F, row.names = F)

##### Combine split datasets #####

split_folder = "./data/rrbs_mouse/thompson/"

# Single
files = dir(split_folder, pattern = "single[^\\.].*", full.names=T)
parts = list()
for (i in 1:length(files)) {
  parts[[i]] = readRDS(files[i])
}
single = Reduce(function(x, y) merge(x, y, by=c(
  "chr", "start", "end"), all=T), parts)
saveRDS(single, paste0(split_folder, "/single.rds"))

# Individual
files = dir(split_folder, pattern = "individual[^\\.].*", full.names=T)
parts = list()
for (i in 1:length(files)) {
  parts[[i]] = readRDS(files[i])
}
individual = Reduce(function(x, y) merge(x, y, by=c(
  "chr", "start", "end", "rep_id", "fragment_id", "strand"), all=T), parts)
saveRDS(individual, paste0(split_folder, "/individual.rds"))

# Collapsed
files = dir(split_folder, pattern = "collapsed[^\\.].*", full.names=T)
parts = list()
for (i in 1:length(files)) {
  parts[[i]] = readRDS(files[i])
}
collapsed = Reduce(function(x, y) merge(x, y, by="clust", all=T), parts)
saveRDS(collapsed, paste0(split_folder, "/collapsed.rds"))

#### COMBINE DATASETS ####

##### Reload collapsed datasets #####

datasets = list()

datasets$stubbs = readRDS("./data/rrbs_mouse/stubbs/collapsed.rds")
datasets$mat_nutrition = readRDS("./data/rrbs_mouse/mat_nutrition/collapsed.rds")
datasets$huntington = readRDS("./data/rrbs_mouse/huntington/collapsed.rds")
datasets$epi_memory = readRDS("./data/rrbs_mouse/epi_memory/collapsed.rds")
datasets$petkovich = readRDS("./data/rrbs_mouse/petkovich/collapsed.rds")
datasets$thompson = readRDS("./data/rrbs_mouse/thompson/collapsed.rds")

meta = list()
# meta$stubbs = read.csv(paste0(paths$out, "/stubbs/meta.txt")) %>%
#   dplyr::select(Run, Sample.Name, Assay.Type, Age, TISSUE, source_name)
tmp = colnames(datasets$stubbs)[seq(2, ncol(datasets$stubbs), 2)]
tmp = str_replace_all(tmp, "numCs_", "")
meta$stubbs = as.data.frame(str_split_fixed(tmp, "_", 3))
colnames(meta$stubbs) = c("MouseID", "Age", "TISSUE")
meta$stubbs$Sample.Name = tmp

meta$mat_nutrition = read.csv("./data/rrbs_mouse/mat_nutrition/meta.txt") %>%
  dplyr::select(Run, Sample.Name, Assay.Type, Age, sex, TISSUE, maternal_diet, adult_diet, Strain)

meta$huntington = read.csv("./data/rrbs_mouse/huntington/meta.txt") %>%
  dplyr::select(Run, Sample.Name, Assay.Type, TISSUE, Genotype) %>%
  mutate(Age = 6*30.5/7) # 6 months old

meta$epi_memory = read.csv("./data/rrbs_mouse/epi_memory/meta.txt") %>%
  dplyr::select(Run, Sample.Name, Assay.Type, Age, gender, TISSUE, Strain, Condition)

meta$petkovich = read.table("./data/rrbs_mouse/petkovich/meta.txt", sep="\t", header=T)
meta$petkovich$Sample_name = toupper(meta$petkovich$Sample_name)
rownames(meta$petkovich) = meta$petkovich$Sample_name

meta$thompson = read.csv("./data/rrbs_mouse/thompson/meta.txt") %>%
  dplyr::select(Sample.Name, Age, sex, Tissue, Strain) %>%
  distinct()

##### Reformat using lists #####

datasets2 = list()
for (ds in names(datasets)) {
  datasets2[[ds]] = list()
  tmp = datasets[[ds]] %>%
    dplyr::filter(!is.na(clust)) %>%
    column_to_rownames("clust")
  nT_inds = grepl("numTs", colnames(tmp))
  nC_inds = grepl("numCs", colnames(tmp))
  datasets2[[ds]]$nTs = tmp[, nT_inds]
  colnames(datasets2[[ds]]$nTs) = str_replace_all(colnames(datasets2[[ds]]$nTs), "numTs_", "")
  datasets2[[ds]]$nCs = tmp[, nC_inds]
  colnames(datasets2[[ds]]$nCs) = str_replace_all(colnames(datasets2[[ds]]$nCs), "numCs_", "")
  if (any(colnames(datasets2[[ds]]$nCs) != colnames(datasets2[[ds]]$nTs))) {
    print("!! Columns not synced between nTs and nCs !!")
  }
}

colnames(datasets2$petkovich$nTs) = toupper(colnames(datasets2$petkovich$nTs))
colnames(datasets2$petkovich$nCs) = toupper(colnames(datasets2$petkovich$nCs))

##### Collapse multiple runs #####

datasets3 = list()
for (ds in names(datasets2)) {
  if (!ds %in% c("mat_nutrition")) {
    datasets3[[ds]] = datasets2[[ds]]
    next
  }
  lookup = meta[[ds]]$Sample.Name
  names(lookup) =  meta[[ds]]$Run
  datasets3[[ds]] = list()
  datasets3[[ds]]$re_info = datasets2[[ds]]$re_info
  nTs = list()
  nCs = list()
  for (subject in unique(lookup)) {
    runs = names(which(lookup == subject))
    run_inds = grepl(paste(runs, collapse="|"), colnames(datasets2$mat_nutrition$nTs))
    if (sum(run_inds) == 0) {
      next
    }
    nT = datasets2$mat_nutrition$nTs[, run_inds]
    nC = datasets2$mat_nutrition$nCs[, run_inds]
    nTs[[subject]] = rowSums(nT)
    nCs[[subject]] = rowSums(nC)
  }
  datasets3[[ds]]$nTs = as.data.frame(nTs)
  datasets3[[ds]]$nCs = as.data.frame(nCs)
}

##### Unite #####

cov = list()
met = list()
for (ds in names(datasets3)) {
  cov[[ds]] = (datasets3[[ds]]$nTs+datasets3[[ds]]$nCs) %>%
    rownames_to_column("clust")
  met[[ds]] = datasets3[[ds]]$nCs %>%
    rownames_to_column("clust")
}

sapply(cov, dim)
cov = Reduce(function(x, y) merge(x, y, by="clust"), cov)
met = Reduce(function(x, y) merge(x, y, by="clust"), met)

# NA means 0
cov[is.na(cov)] = 0
met[is.na(met)] = 0

##### Remove low depth samples #####

depth = colSums(cov[,-c(1)])
hist(log10(depth), breaks=100)
min_depth = 1e6
abline(v=log10(min_depth), col="red")

cov = cov[, c(T, depth > min_depth)]
met = met[, c(T, depth > min_depth)]

##### Remove low coverage REs #####

th = 100
prev = 0.9

met = column_to_rownames(met, "clust")
cov = column_to_rownames(cov, "clust")

passing = rowSums(cov > th) > prev * ncol(cov)
table(passing)
met2 = met[passing,]
cov2 = cov[passing,]

data_rrbs = list(
  met = met,
  cov = cov,
  met_hc = met2,
  cov_hc = cov2
)

##### Make re_info #####

re_order = c("LINE", "SINE", "LTR", "DNA", "DNA?", "Retroposon", "RC", 
             "tRNA", "srpRNA", "snRNA", "scRNA", "rRNA", 
             "Satellite", "Simple_repeat", "Low_complexity",
             "Unspecified", "Unknown")
fam_ids = read.csv(paths$re_ids, row.names = 1) %>%
  rownames_to_column("fam_id")%>%
  mutate(class = factor(class, levels=re_order))%>%
  mutate(superf = factor(superf, levels = unique(.$superf[base::order(.$class, .$superf)]))) %>%
  dplyr::arrange(class)
re_info = str_split(rownames(met), pattern="_", simplify = T)
re_info = data.frame(re_info) %>%
  mutate(feature = rownames(met)) %>%
  dplyr::rename(fam_id="X1", annotation3="X2") %>%
  merge(., fam_ids, by="fam_id") %>%
  column_to_rownames("feature")
data_rrbs$re_info = re_info[rownames(cov), ]

##### Unite meta #####

meta_rrbs = data.frame(
  study = character(),
  individual = character(),
  tissue = character(),
  age = numeric(),
  sex = character(),
  strain = character(),
  treatment = character()
)
meta_rrbs = meta$stubbs %>%
  column_to_rownames(var = "Sample.Name") %>%
  mutate(age = as.numeric(str_replace_all(Age, "wk", ""))) %>%
  mutate(tissue = tolower(TISSUE)) %>%
  mutate(individual = MouseID) %>%
  mutate(sex = "M") %>%
  mutate(strain = "C57BL/6-BABR") %>%
  mutate(study = "stubbs") %>%
  mutate(treatment = "untreated") %>%
  dplyr::select(study, individual, tissue, age, sex, strain, treatment) %>%
  rbind(meta_rrbs, .)
meta_rrbs = meta$mat_nutrition %>%
  dplyr::select(-Run) %>%
  distinct() %>%
  column_to_rownames(var = "Sample.Name") %>%
  mutate(age = 9) %>%
  mutate(tissue = TISSUE) %>%
  mutate(individual = rownames(.)) %>%
  mutate(sex = fct_recode(sex, "M" = "male", "F" = "female")) %>%
  mutate(strain = Strain) %>%
  mutate(study = "mat_nutrition") %>%
  mutate(treatment = paste0("mat ", maternal_diet, "-adult ", adult_diet)) %>%
  dplyr::select(study, individual, tissue, age, sex, strain, treatment) %>%
  rbind(meta_rrbs, .)
meta_rrbs = meta$huntington %>%
  column_to_rownames(var = "Sample.Name")%>%
  mutate(study = "huntington") %>%
  mutate(individual = rownames(.)) %>%
  mutate(tissue = tolower(TISSUE)) %>%
  mutate(age = Age) %>%
  mutate(sex = NA) %>%
  mutate(strain="huntington model") %>%
  mutate(treatment = Genotype) %>%
  dplyr::select(study, individual, tissue, age, sex, strain, treatment) %>%
  rbind(meta_rrbs, .)
meta_rrbs = meta$epi_memory %>%
  dplyr::filter(Assay.Type == "Bisulfite-Seq") %>%
  dplyr::filter(Strain == "C57/bl") %>%
  column_to_rownames(var="Run") %>%
  mutate(study = "epi_memory") %>%
  mutate(individual = Sample.Name) %>%
  mutate(tissue = tolower(TISSUE)) %>%
  mutate(age = as.numeric(str_replace_all(Age, "w", ""))) %>%
  mutate(sex = fct_recode(gender, "M" = "male", "F" = "female")) %>%
  mutate(strain="C57BL/6") %>%
  mutate(treatment = fct_recode(Condition, "untreated" = "NORMAL")) %>%
  dplyr::select(study, individual, tissue, age, sex, strain, treatment) %>%
  rbind(meta_rrbs, .)
meta_rrbs = meta$petkovich %>%
  mutate(study = "petkovich") %>%
  mutate(individual = Sample_name) %>%
  mutate(tissue = tolower(Tissue)) %>%
  mutate(age = 4.3 * Age_m) %>%
  mutate(sex = Sex) %>%
  mutate(strain = Strain) %>%
  mutate(treatment = paste(Genotype, Diet, sep="-"))%>%
  dplyr::select(study, individual, tissue, age, sex, strain, treatment) %>%
  rbind(meta_rrbs, .)
meta_rrbs = meta$thompson %>%
  dplyr::filter(Age != "") %>%
  mutate(study = "thompson") %>%
  mutate(individual = Sample.Name) %>%
  mutate(tissue = tolower(Tissue)) %>%
  mutate(age = 4.3 * as.numeric(gsub("mo", "", Age))) %>%
  mutate(sex = fct_recode(sex, "M" = "male", "F" = "female")) %>%
  mutate(strain = Strain) %>%
  mutate(treatment = "untreated") %>%
  column_to_rownames(var="Sample.Name") %>%
  dplyr::select(study, individual, tissue, age, sex, strain, treatment) %>%
  rbind(meta_rrbs, .)
    
csamps = intersect(rownames(meta_rrbs), colnames(data_rrbs$met))
setdiff(rownames(meta_rrbs), colnames(data_rrbs$met))
setdiff(colnames(data_rrbs$met), rownames(meta_rrbs))
data_rrbs$met = data_rrbs$met[, csamps]
data_rrbs$cov = data_rrbs$cov[, csamps]
data_rrbs$met_hc = data_rrbs$met_hc[, csamps]
data_rrbs$cov_hc = data_rrbs$cov_hc[, csamps]
meta_rrbs = meta_rrbs[csamps, ]
rownames(meta_rrbs) = paste(meta_rrbs$study, rownames(meta_rrbs), sep=".")
colnames(data_rrbs$met) = rownames(meta_rrbs)
colnames(data_rrbs$cov) = rownames(meta_rrbs)
colnames(data_rrbs$met_hc) = rownames(meta_rrbs)
colnames(data_rrbs$cov_hc) = rownames(meta_rrbs)

# <Output>
save(data_rrbs, meta_rrbs, file=paste0(paths$out, "/rrbs_data.Rdata"))
