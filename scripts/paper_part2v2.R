library(tidyverse)
library(ggpubr)
library(limma)
library(patchwork)
library(caret)
library(glmnet)
library(ModelMetrics)
library(ggdendro)
library(clevr)
library(rasterpdf)
library(data.table)
library(survival)
library(survminer)
library(metafor)
library(survminer)
library(randomForestSRC)
library(nestedcv)
library(caret)
library(pracma)
library(gtools)

setwd("/scratch/fmorandi/internal/RE_clock/")

source("./scripts/utils.R")

paths = list()
paths$data = "./data"
paths$out = "./results/final"

#### PLOTTING SETTINGS ####

w = 210 - 25.4*2 # was 174
h = 297 - 25.4*2 # was 230
w_in = w*0.0393701
h_in = h*0.0393701

theme_set(theme_light(base_size = 8))
update_geom_defaults("text", list(size = 8*0.35))

re_order = c("LINE", "SINE", "LTR", "DNA", "DNA?", "Retroposon", "RC", 
             "tRNA", "srpRNA", "snRNA", "scRNA", "rRNA", 
             "Satellite", "Simple_repeat", "Low_complexity",
             "Unspecified", "Unknown")

#### FUNCTIONS ####

simple_split2 = function(data, meta, re_info, alphas) {
  clock = simple_split(as.matrix(data), ageT(meta$age, 20),
                       cv_function = cv_glmnet, cv_args=list(alphas=alphas),
                       filt_function=filt_by_min_n, filt_args = list(minn=5, re_info=re_info))
  return(clock)
}
predict_plot2 = function(clock, meta, title) {
  p=predict_plot(clock$preds, meta, "study", untransform_age = 20)+
    ggtitle(title, subtitle = sprintf("Lambda = %.3f, Alpha = %.3f", clock$fit$lambda, clock$fit$alpha))
  return(p)
}

#### LOAD DATA ####

# Load tables
data_infi = readRDS(paste0(paths$data, "/processed/infinium_data.rds"))
load(paste0(paths$data, "/processed/infinium_data_collapsed1.Rdata"))
load(paste0(paths$data, "/processed/infinium_data_collapsed2.Rdata"))
load(paste0(paths$data, "/processed/infinium_data_collapsed3.Rdata"))

# Family IDs for collapsing WHI data
fam_ids = read.csv(paste0(paths$data, "/annotation/GRCh37_repIds.csv"), row.names = 1) %>%
  rownames_to_column("fam_id")

# Probe info
probe_info = readRDS(paste0(paths$data, "/processed/infinium_probes.rds"))
rownames(probe_info) = probe_info$IlmnID
probe_info = probe_info[colnames(data_infi), ]

# Metadata
meta_infi = fread(paste0(paths$data, "/processed/infinium_meta.tsv"), data.table=F, sep="\t")
meta_infi = meta_infi[, -c(1)]
rownames(meta_infi) = meta_infi$ID

# Collapsed repeat info
re_info1 = re_info1 %>%
  mutate(class = factor(class, levels=re_order))%>%
  mutate(superf = factor(superf, levels = unique(.$superf[base::order(.$class, .$superf)]))) %>%
  dplyr::arrange(class)
re_info2 = re_info2 %>%
  mutate(class = factor(class, levels=re_order))%>%
  mutate(superf = factor(superf, levels = unique(.$superf[base::order(.$class, .$superf)]))) %>%
  dplyr::arrange(class)
re_info3 = re_info3 %>%
  mutate(class = factor(class, levels=re_order))%>%
  mutate(superf = factor(superf, levels = unique(.$superf[base::order(.$class, .$superf)]))) %>%
  dplyr::arrange(class)

# Published clocks
clocks = list()
clocks$hannum = read.table(paste0(paths$data, "/published_clocks/hannum_blood_coefs.txt"), header=T)
rownames(clocks$hannum) = clocks$hannum$IlmnID
clocks$horvath_pt = read.table(paste0(paths$data, "/published_clocks/horvath_pantissue_coefs.txt"), header=T)
rownames(clocks$horvath_pt) = clocks$horvath_pt$IlmnID
clocks$horvath_sb = read.table(paste0(paths$data, "/published_clocks/horvath_skinblood_coefs.txt"), header=T)
rownames(clocks$horvath_sb) = clocks$horvath_sb$IlmnID
clocks$horvath_pheno = read.table(paste0(paths$data, "/published_clocks/horvath_pheno_coefs.txt"), header=T)
rownames(clocks$horvath_pheno) = clocks$horvath_pheno$IlmnID

#### FIGURE STORAGE ####

load(paste0(paths$out, "/figures.Rdata"))

#### OUTLIER REMOVAL ####

outlier1 = outlier_detection(data_re1, meta_infi, subgroups = c("sex", "health"),
                             savepath = paste0(paths$out, "/outliers1.pdf"))
outlier2 = outlier_detection(data_re2, meta_infi, subgroups = c("sex", "health"),
                             savepath = paste0(paths$out, "/outliers2.pdf"))
outlier3 = outlier_detection(data_re3, meta_infi, subgroups = c("sex", "health"),
                             savepath = paste0(paths$out, "/outliers3.pdf"))
meta_infi$outlier = outlier1 | outlier2 | outlier3
meta_infi = subset(meta_infi, !meta_infi$outlier)
data_re1 = data_re1[rownames(meta_infi), ]
data_re2 = data_re2[rownames(meta_infi), ]
data_re3 = data_re3[rownames(meta_infi), ]
data_infi = data_infi[rownames(meta_infi), ]

#### RESTRICT TO SELFISH ####

# Define selfish
probe_info$selfish = ifelse(probe_info$class %in% c("LINE", "SINE", "LTR", "DNA"), "Selfish", "Non-selfish")
re_info1$selfish = ifelse(re_info1$class %in% c("LINE", "SINE", "LTR", "DNA"), "Selfish", "Non-selfish")
re_info2$selfish = ifelse(re_info2$class %in% c("LINE", "SINE", "LTR", "DNA"), "Selfish", "Non-selfish")
re_info3$selfish = ifelse(re_info3$class %in% c("LINE", "SINE", "LTR", "DNA"), "Selfish", "Non-selfish")

# Restrict
cpg_selfish = probe_info[probe_info$selfish == "Selfish", "IlmnID"]
# cpg_selfish_nopro = intersect(cpg_selfish, probe_info[!probe_info$annotation2 %in% c("Promoter", "Exon"), "IlmnID"])
re1_selfish = rownames(re_info1)[re_info1$selfish == "Selfish"]
re2_selfish = rownames(re_info2)[re_info2$selfish == "Selfish"]
re3_selfish = str_replace(rownames(re_info3)[re_info3$selfish == "Selfish"], " ", ".")
data_selfish = data_infi[, cpg_selfish]
# data_selfish_nopro = data_selfish[, cpg_selfish_nopro]
data_re1_selfish = data_re1[, re1_selfish]
data_re2_selfish = data_re2[, re2_selfish]
data_re3_selfish = data_re3[, re3_selfish]

# Feature numbers
table(probe_info$selfish)
table(re_info1$selfish)
table(re_info2$selfish)
table(re_info3$selfish)
table(re_info1[re_info1$n >= 5, "selfish"])
table(re_info2[re_info2$n >= 5, "selfish"])
table(re_info3[re_info3$n >= 5, "selfish"])

# Check
all.equal(rownames(data_selfish), rownames(meta_infi))
all.equal(rownames(data_re1_selfish), rownames(meta_infi))
all.equal(rownames(data_re2_selfish), rownames(meta_infi))
all.equal(rownames(data_re3_selfish), rownames(meta_infi))

#### CLOCK TRAINING ####

##### Settings #####

# Universal settings
healthy_ind = which(meta_infi$health == "Healthy")
not_GSE64495_ind = intersect(which(meta_infi$study != "GSE64495"), healthy_ind)
# alphas = c(0, 0.4, 0.5, 0.8, 0.9, 0.95, 1)
alphas = seq(0, 1, 0.1)

clocks = readRDS(paste0(paths$out, "/clocks.rds"))

##### Single CpG, selfish only #####

# === Single CpG, selfish, any context ===
# clocks$single_cpg = simple_split(as.matrix(data_selfish[not_GSE64495_ind, ]), ageT(meta_infi[not_GSE64495_ind, "age"], 20),
#                                           cv_function = cv_glmnet, cv_args=list(alphas=alphas))
fig$M3$clock_single_cpg = predict_plot(clocks$single_cpg$preds, meta_infi, "study", untransform_age = 20)+
  ggtitle("Single CpG clock, excludes GSE64495",
          subtitle = sprintf("Lambda = %.3f, Alpha = %.3f",
                             clocks$single_cpg$fit$lambda,
                             clocks$single_cpg$fit$alpha))
ggsave(plot=fig$M3$clock_single_cpg, paste0(paths$out, "/clock_single_cpg.png"), width=0.6*w, height=0.3*h, units="mm", dpi = 300)

# # === Single CpG, selfish, exclude promoters and exons ===
# clocks$single_cpg_nopro = simple_split(as.matrix(data_selfish_nopro[not_GSE64495_ind, ]), ageT(meta_infi[not_GSE64495_ind, "age"], 20),
#                                       cv_function = cv_glmnet, cv_args=list(alphas=alphas))
# predict_plot(clocks$single_cpg_nopro$preds, meta_infi, "study", untransform_age = 20)+
#   ggtitle("Single CpG clock, no promoters, excludes GSE64495",
#           subtitle = sprintf("Lambda = %.3f, Alpha = %.3f",
#                              clocks$single_cpg_nopro$fit$lambda,
#                              clocks$single_cpg_nopro$fit$alpha))

##### Combined, selfish only #####

# === Combined (no genic context), min 5 ===
# clocks$combined1_min5 = simple_split2(data_re1_selfish[not_GSE64495_ind, ], meta_infi[not_GSE64495_ind, ], re_info1[re1_selfish, ], alphas)
predict_plot2(clocks$combined1_min5, meta_infi, "Re1 clock, minn=5, excludes GSE64495, selfish only")
ggsave(paste0(paths$out, "/clock_combined1_min5.png"), width=0.6*w, height=0.3*h, units="mm", dpi = 300)

# === Combined (broad genic context), min 5 ===
# clocks$combined2_min5 = simple_split2(data_re2_selfish[not_GSE64495_ind, ], meta_infi[not_GSE64495_ind, ], re_info2[re2_selfish, ], alphas)
fig$M3$clock_combined = predict_plot2(clocks$combined2_min5, meta_infi, "Re2 clock, minn=5, excludes GSE64495, selfish only")
ggsave(plot=fig$M3$clock_combined, paste0(paths$out, "/clock_combined2_min5.png"), width=0.6*w, height=0.3*h, units="mm", dpi = 300)

# === Combined (detailed genic context), min 5 ===
# clocks$combined3_min5 = simple_split2(data_re3_selfish[not_GSE64495_ind, ], meta_infi[not_GSE64495_ind, ], re_info3[re3_selfish, ], alphas)
predict_plot2(clocks$combined3_min5, meta_infi, "Re3 clock, minn=5, excludes GSE64495, selfish only")
ggsave(paste0(paths$out, "/clock_combined3_min5.png"), width=0.6*w, height=0.3*h, units="mm", dpi = 300)

##### Save #####

# saveRDS(clocks, file = paste0(paths$out, "/clocks.rds"))

#### CLOCK FEATURES ####

# Get feature annotation
feats_single = data.frame(as.matrix(coef(clocks$single_cpg$fit))) %>%
  dplyr::filter(s0 != 0) %>%
  rownames_to_column("IlmnID") %>%
  merge(probe_info, ., by="IlmnID", all.y=T)
# feats_single_nopro = data.frame(as.matrix(coef(clocks$single_cpg_nopro$fit))) %>%
#   dplyr::filter(s0 != 0) %>%
#   rownames_to_column("IlmnID") %>%
#   merge(probe_info, ., by="IlmnID", all.y=T)
feats_aggreg = data.frame(as.matrix(coef(clocks$combined2_min5$fit))) %>%
  dplyr::filter(s0 != 0) %>%
  rownames_to_column("reID") %>%
  merge(re_info2, ., by.x=0, by.y="reID", all.y=T) %>%
  dplyr::rename("reID" = "Row.names")

# Save (Plot made together with RRBS clock)
write.csv(feats_single, file=paste0(paths$out, "/clock_feats_single.csv"))
# write.csv(feats_single_nopro, file=paste0(paths$out, "/clock_feats_single_nopro.csv"))
write.csv(feats_aggreg, file=paste0(paths$out, "/clock_feats_combined.csv"))

# Summary by TE class
tmp1 = feats_single %>%
  dplyr::filter(!is.na(class)) %>%
  group_by(class) %>%
  dplyr::summarize(nfeat=n()) %>%
  mutate(clock = "Individual CpG\nTE clock")
write.csv(tmp1, paste0(paths$out, "/class_composition_single.csv"))
# tmp1 = feats_single_nopro %>%
#   dplyr::filter(!is.na(class)) %>%
#   group_by(class) %>%
#   dplyr::summarize(nfeat=n()) %>%
#   mutate(clock = "Individual CpG\nTE clock")
# write.csv(tmp1, paste0(paths$out, "/class_composition_single_nopro.csv"))
tmp2 = feats_aggreg %>%
  dplyr::filter(!is.na(class)) %>%
  group_by(class) %>%
  dplyr::summarize(nfeat=n()) %>%
  mutate(clock = "Combined CpG\nTE clock")
write.csv(tmp2, paste0(paths$out, "/class_composition_aggreg.csv"))

# Intersection with published clocks
intersect(feats_single$IlmnID, clocks$hannum$IlmnID)
intersect(feats_single$IlmnID, clocks$horvath_pt$IlmnID)
intersect(feats_single$IlmnID, clocks$horvath_sb$IlmnID)
intersect(feats_single$IlmnID, clocks$horvath_pheno$IlmnID)
#
# intersect(feats_single_nopro$IlmnID, clocks$hannum$IlmnID)
# intersect(feats_single_nopro$IlmnID, clocks$horvath_pt$IlmnID)
# intersect(feats_single_nopro$IlmnID, clocks$horvath_sb$IlmnID)
# intersect(feats_single_nopro$IlmnID, clocks$horvath_pheno$IlmnID)

clocks_oi = c("hannum", "horvath_pt", "horvath_sb", "horvath_pheno")
res_intersection = data.frame()
# Intersections between published clocks
for (i in 1:(length(clocks_oi)-1)) {
  for (j in (i+1):length(clocks_oi)) {
    res_intersection[nrow(res_intersection)+1, "clock1"] = clocks_oi[i]
    res_intersection[nrow(res_intersection), "clock2"] = clocks_oi[j]
    res_intersection[nrow(res_intersection), "N1"] = nrow(clocks[[clocks_oi[i]]])-1
    res_intersection[nrow(res_intersection), "N2"] = nrow(clocks[[clocks_oi[j]]])-1
    common = intersect(clocks[[clocks_oi[i]]]$IlmnID, clocks[[clocks_oi[j]]]$IlmnID)
    common = setdiff(common, "(Intercept)")
    res_intersection[nrow(res_intersection), "Ncommon"] = length(common)
    res_intersection[nrow(res_intersection), "common"] = paste(common, collapse=", ")
  }
}
# res_intersection_nopro = res_intersection
# Intersection between published clocks and single cpg TE clock
for (i in 1:(length(clocks_oi))) {
    res_intersection[nrow(res_intersection)+1, "clock1"] = "single_cpg"
    res_intersection[nrow(res_intersection), "clock2"] = clocks_oi[i]
    res_intersection[nrow(res_intersection), "N1"] = nrow(feats_single)-1
    res_intersection[nrow(res_intersection), "N2"] = nrow(clocks[[clocks_oi[i]]])-1
    common = intersect(feats_single$IlmnID, clocks[[clocks_oi[i]]]$IlmnID)
    common = setdiff(common, "(Intercept)")
    res_intersection[nrow(res_intersection), "Ncommon"] = length(common)
    res_intersection[nrow(res_intersection), "common"] = paste(common, collapse=", ")
}
# # Intersection between published clocks and single cpg TE clock with no promoters
# for (i in 1:(length(clocks_oi))) {
#   res_intersection_nopro[nrow(res_intersection_nopro)+1, "clock1"] = "single_cpg_nopro"
#   res_intersection_nopro[nrow(res_intersection_nopro), "clock2"] = clocks_oi[i]
#   res_intersection_nopro[nrow(res_intersection_nopro), "N1"] = nrow(feats_single_nopro)-1
#   res_intersection_nopro[nrow(res_intersection_nopro), "N2"] = nrow(clocks[[clocks_oi[i]]])-1
#   common = intersect(feats_single_nopro$IlmnID, clocks[[clocks_oi[i]]]$IlmnID)
#   common = setdiff(common, "(Intercept)")
#   res_intersection_nopro[nrow(res_intersection_nopro), "Ncommon"] = length(common)
#   res_intersection_nopro[nrow(res_intersection_nopro), "common"] = paste(common, collapse=", ")
# }
# Extra stats
res_intersection = res_intersection %>%
  mutate(Nunion = N1 + N2 - Ncommon) %>%
  mutate(IoU = Ncommon/Nunion)
# res_intersection_nopro = res_intersection_nopro %>%
#   mutate(Nunion = N1 + N2 - Ncommon) %>%
#   mutate(IoU = Ncommon/Nunion)

# Plot single CpG
feats_single = feats_single %>%
  mutate(shared_hannum = IlmnID %in% clocks$hannum$IlmnID) %>%
  mutate(shared_horvath_pt = IlmnID %in% clocks$horvath_pt$IlmnID) %>%
  mutate(shared_horvath_sb = IlmnID %in% clocks$horvath_sb$IlmnID) %>%
  mutate(shared_horvath_pheno = IlmnID %in% clocks$horvath_pheno$IlmnID)
tmp = feats_single %>%
  pivot_longer(cols=contains("shared"), names_to = "comp", values_to = "shared") %>%
  dplyr::filter(shared) %>%
  dplyr::filter(IlmnID != "(Intercept)")
xrange = range(feats_single$s0)*1.1
p1 = ggplot(feats_single, aes(s0))+
  geom_histogram(bins=100)+
  scale_y_sqrt()+
  xlim(xrange)
p2 = ggplot(tmp, aes(x=s0, y=comp, group=IlmnID))+
  geom_point()+
  geom_line()+
  xlim(xrange)
fig$S5$clock_overlap = p1/p2&theme(axis.title.y=element_blank())&labs(x="Coefficient")
ggsave(plot = fig$S5$clock_overlap, paste0(paths$out, "/clock_overlap.png"), width=1*w, height=0.5*h, units="mm", dpi = 300)

# # Plot single CpG no promoter
# feats_single_nopro = feats_single_nopro %>%
#   mutate(shared_hannum = IlmnID %in% clocks$hannum$IlmnID) %>%
#   mutate(shared_horvath_pt = IlmnID %in% clocks$horvath_pt$IlmnID) %>%
#   mutate(shared_horvath_sb = IlmnID %in% clocks$horvath_sb$IlmnID) %>%
#   mutate(shared_horvath_pheno = IlmnID %in% clocks$horvath_pheno$IlmnID)
# tmp = feats_single_nopro %>%
#   pivot_longer(cols=contains("shared"), names_to = "comp", values_to = "shared") %>%
#   dplyr::filter(shared) %>%
#   dplyr::filter(IlmnID != "(Intercept)")
# xrange = range(feats_single_nopro$s0)*1.1
# p1 = ggplot(feats_single_nopro, aes(s0))+
#   geom_histogram(bins=100)+
#   scale_y_sqrt()+
#   xlim(xrange)
# p2 = ggplot(tmp, aes(x=s0, y=comp, group=IlmnID))+
#   geom_point()+
#   geom_line()+
#   xlim(xrange)
# fig$S3$clock_overlap_nopro = p1/p2&theme(axis.title.y=element_blank())&labs(x="Coefficient")
# ggsave(plot = fig$S3$clock_overlap_nopro, paste0(paths$out, "/clock_overlap_nopro.png"), width=1*w, height=0.5*h, units="mm", dpi = 300)

#### CLOCK COMPARISON ####

# datasets = unique(meta_infi$study)
datasets = "GSE64495"
cpgs_oi = union(clocks$hannum$IlmnID, clocks$horvath_pt$IlmnID)
cpgs_oi = union(cpgs_oi, clocks$horvath_sb$IlmnID)
cpgs_oi = union(cpgs_oi, clocks$horvath_pheno$IlmnID)
cpgs_oi = setdiff(cpgs_oi, "(Intercept)")
clocks_oi = c("hannum", "horvath_pt", "horvath_sb")
# my_clocks = c("single_not_GSE64495", "re2_minn5_not_GSE64495")

res = list(
  "rmse" = c(),
  "mae" = c(),
  "r" = c(),
  "test_dataset" = c(),
  "clock_used" = c()
)

for (ds in datasets) {
  load(paste0(paths$data, "/infinium/", ds, "/", ds, ".Rdata"))
  if(!"ID_REF" %in% colnames(data)) {
    data$ID_REF = rownames(data)
  }
  data = as.data.frame(subset(data, ID_REF %in% cpgs_oi))
  rownames(data) = data$ID_REF
  data$ID_REF = NULL
  data = data %>%
    mutate_if(is.character, as.numeric) %>%
    t()
  data = data[intersect(rownames(data), rownames(data_re1)), ]
  labels = meta_infi[rownames(data), "age"]
  # === Test published single cpg clocks ===
  for (cl in clocks_oi) {
    if (length(setdiff(clocks[[cl]]$IlmnID, colnames(data))) > 1) {
      print(sprintf("Missing features, skipping %s for clock %s", ds, cl))
      next
    }
    untr = cl %in% c("horvath_pt", "horvath_sb")
    their_preds = predict_from_coefs(clocks[[cl]], data, untransform_age = untr)
    res$rmse = c(res$rmse, rmse(labels, their_preds))
    res$mae = c(res$mae, mae(labels, their_preds))
    res$r = c(res$r, cor(labels, their_preds))
    res$test_dataset = c(res$test_dataset, ds)
    res$clock_used = c(res$clock_used, cl)
  }
  # === Test my clocks ===
  # Single CpG
  my_coefs = feats_single %>%
    dplyr::rename("Coef" = "s0")
  rownames(my_coefs) = my_coefs$IlmnID
  my_preds = predict_from_coefs(my_coefs, data_selfish[meta_infi$study == ds, ], untransform_age = 20)
  res$rmse = c(res$rmse, rmse(labels, my_preds))
  res$mae = c(res$mae, mae(labels, my_preds))
  res$r = c(res$r, cor(labels, my_preds))
  res$test_dataset = c(res$test_dataset, ds)
  res$clock_used = c(res$clock_used, "single_not_GSE64495")
  # Single CpG no promoters
  # my_coefs = feats_single_nopro %>%
  #   dplyr::rename("Coef" = "s0")
  # rownames(my_coefs) = my_coefs$IlmnID
  # my_preds = predict_from_coefs(my_coefs, data_selfish_nopro[meta_infi$study == ds, ], untransform_age = 20)
  # res$rmse = c(res$rmse, rmse(labels, my_preds))
  # res$mae = c(res$mae, mae(labels, my_preds))
  # res$r = c(res$r, cor(labels, my_preds))
  # res$test_dataset = c(res$test_dataset, ds)
  # res$clock_used = c(res$clock_used, "single_nopro")
  # Aggregated
  my_coefs = feats_aggreg %>%
    dplyr::rename("Coef" = "s0")
  rownames(my_coefs) = my_coefs$reID
  my_preds = predict_from_coefs(my_coefs, data_re2_selfish[meta_infi$study == ds, ], untransform_age = 20)
  res$rmse = c(res$rmse, rmse(labels, my_preds))
  res$mae = c(res$mae, mae(labels, my_preds))
  res$r = c(res$r, cor(labels, my_preds))
  res$test_dataset = c(res$test_dataset, ds)
  res$clock_used = c(res$clock_used, "re2_min5_not_GSE64495")
}

res = as.data.frame(res)
res$clock_used = factor(res$clock_used, levels=res$clock_used)
p1 = ggplot(res, aes(x=clock_used, y=rmse, fill=clock_used))+
  geom_bar(stat="identity")+
  scale_y_continuous(expand=expansion(mult=c(0, 0.1)))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  ggtitle("RMSE")
p2 = ggplot(res, aes(x=clock_used, y=mae, fill=clock_used))+
  geom_bar(stat="identity")+
  scale_y_continuous(expand=expansion(mult=c(0, 0.1)))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  ggtitle("MAE")
p3 = ggplot(res, aes(x=clock_used, y=r, fill=clock_used))+
  geom_bar(stat="identity")+
  scale_y_continuous(expand=expansion(mult=c(0, 0.1)))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  ggtitle("r")
fig$M3$benchmark_chrono = p1+p2+p3+plot_layout(guides="collect")+plot_annotation(title="Test on GSE64495")
ggsave(plot=fig$M3$benchmark_chrono, paste0(paths$out, "/comparison_on_GSE64495.png"), width=1*w, height=0.3*h, units="mm", dpi = 300)

#### CHECKPOINT ####

save(fig, file=paste0(paths$out, "/figures.Rdata"))
# save.image(paste0(paths$out, "/checkpoint_p2.Rdata"))
load(paste0(paths$out, "/checkpoint_p2.Rdata"))

#### WHI DATA ####

##### Load ####

load(paste0(paths$data, "/processed/whi_data.Rdata"))
meta_whi = read.table(paste0(paths$data, "/processed/whi_meta.tsv"))

rownames(data_whi) = data_whi$IlmnID
data_whi$IlmnID = NULL

data_whi = mutate_if(data_whi, is.character, as.numeric)
rownames(meta_whi) = meta_whi$SampleID
meta_whi = subset(meta_whi, !is.na(AGE))
common = intersect(colnames(data_whi), meta_whi$SampleID)
data_whi = data_whi[, common]

##### Collapse ##### 

common_cpgs = intersect(colnames(data_infi), rownames(data_whi))
probe_info_whi = merge(probe_info[common_cpgs, ], fam_ids, by=c("class", "superf", "fam"), all.x=T)
rownames(probe_info_whi) = probe_info_whi$IlmnID
probe_info_whi = probe_info_whi[common_cpgs, ]

re_info_whi1 = probe_info_whi %>%
  mutate(clust = fam_id)
data_re1_whi = collapse_to_clusters(data_whi[common_cpgs, ], re_info_whi1$clust, rows_are_features=T)
re_info_whi1 = re_info_whi1 %>%
  group_by(class, superf, fam, fam_id, clust) %>%
  dplyr::filter(!is.na(clust)) %>%
  dplyr::summarize(n=n())%>%
  column_to_rownames("clust")

re_info_whi2 = probe_info_whi %>%
  mutate(clust = paste(fam_id, annotation3, sep="_"))
re_info_whi2[is.na(re_info_whi2$class), "clust"] = NA
data_re2_whi = collapse_to_clusters(data_whi[common_cpgs, ], re_info_whi2$clust, rows_are_features=T)
re_info_whi2 = re_info_whi2 %>%
  group_by(class, superf, fam, fam_id, annotation3, clust) %>%
  dplyr::filter(!is.na(clust)) %>%
  dplyr::summarize(n=n())%>%
  column_to_rownames("clust")

# re_info_whi3 = probe_info_whi %>%
#   mutate(clust = paste(fam_id, annotation2, sep="_"))
# re_info_whi3[is.na(re_info_whi3$class), "clust"] = NA
# data_re3_whi = collapse_to_clusters(data_whi[common_cpgs, ], re_info_whi3$clust, rows_are_features=T)
# re_info_whi3 = re_info_whi3 %>%
#   group_by(class, superf, fam, fam_id, annotation2, clust) %>%
#   dplyr::filter(!is.na(clust)) %>%
#   dplyr::summarize(n=n())%>%
#   column_to_rownames("clust")

##### Keep CpGs of interest #####

# Published clock cpgs
cpgs_oi = union(clocks$hannum$IlmnID, clocks$horvath_pt$IlmnID)
cpgs_oi = union(cpgs_oi, clocks$horvath_sb$IlmnID)
cpgs_oi = union(cpgs_oi, clocks$horvath_pheno$IlmnID)
cpgs_oi = union(cpgs_oi, feats_single$IlmnID)
cpgs_oi = setdiff(cpgs_oi, "(Intercept)")

probe_info_yl1 = probe_info %>%
  dplyr::filter(superf == "LINE/L1") %>%
  dplyr::filter(fam %in% c("L1HS", "L1PA2", "L1PA3", "L1PA4"))
probe_info_the1 = probe_info %>%
  dplyr::filter(superf == "LTR/ERVL-MaLR") %>%
  dplyr::filter(fam %in% c("THE1A", "THE1C"))
cpgs_oi = union(cpgs_oi, probe_info_yl1$IlmnID)
cpgs_oi = union(cpgs_oi, probe_info_the1$IlmnID)

data_whi = data_whi[rownames(data_whi) %in% cpgs_oi, ]
data_whi = data.frame(t(data_whi))

data_re1_whi = data.frame(t(data_re1_whi))
data_re2_whi = data.frame(t(data_re2_whi))
# data_re3_whi = data.frame(t(data_re3_whi))

##### Checkpoint #####

# save(data_re1_whi, re_info_whi1, data_re2_whi, re_info_whi2, data_whi, meta_whi, file=paste0(paths$out, "/checkpoint_whi.Rdata"))
load(paste0(paths$out, "/checkpoint_whi.Rdata"))

data_re1_whi = data_re1_whi[, colnames(data_re1)]
data_re2_whi = data_re2_whi[, colnames(data_re2)]
# data_re3_whi = data_re3_whi[, colnames(data_re3)]

rownames(data_re1_whi) = rownames(data_whi)
rownames(data_re2_whi) = rownames(data_whi)
# rownames(data_re3_whi) = rownames(data_whi)

data_whi = data_whi[as.character(meta_whi$SampleID), ]
data_re1_whi = data_re1_whi[as.character(meta_whi$SampleID), ]
data_re2_whi = data_re2_whi[as.character(meta_whi$SampleID), ]
# data_re3_whi = data_re3_whi[as.character(meta_whi$SampleID), ]

##### Compare predictions #####

sum(is.na(data_whi)) / (nrow(data_whi) * ncol(data_whi))
data_whi = data.frame(makeX(data_whi, na.impute = T))

res = list(
  "rmse" = c(),
  "mae" = c(),
  "r" = c(),
  "clock_used" = c()
)
clocks_oi = c("hannum", "horvath_pt", "horvath_sb", "horvath_pheno")
for (cl in clocks_oi) {
  tmp = data_whi[, colnames(data_whi) %in% clocks[[cl]]$IlmnID]
  # tmp = tmp[rownames(data_whi)[rowSums(is.na(tmp)) == 0], ]
  if (length(setdiff(clocks[[cl]]$IlmnID, colnames(tmp))) > 1) {
    print(sprintf("Missing features, skipping %s", cl))
    next
  }
  untr = cl %in% c("horvath_pt", "horvath_sb")
  if (any(rownames(tmp) != rownames(meta_whi))) {
    print("tmp rows don't match meta_whi")
  }
  meta_whi[[cl]] = as.numeric(predict_from_coefs(clocks[[cl]], tmp, untransform_age = untr))
  tmp = meta_whi %>%
    drop_na(cl)
  res$rmse = c(res$rmse, rmse(tmp$AGE, tmp[[cl]]))
  res$mae = c(res$mae, mae(tmp$AGE, tmp[[cl]]))
  res$r = c(res$r, cor(tmp$AGE, tmp[[cl]]))
  res$clock_used = c(res$clock_used, cl)
}
# === Test my clocks ===
# Single CpG
my_coefs = feats_single %>%
  dplyr::rename("Coef" = "s0")
rownames(my_coefs) = my_coefs$IlmnID
meta_whi$single = my_preds = predict_from_coefs(my_coefs, data_whi, untransform_age = 20)
tmp = meta_whi %>%
  drop_na(single)
res$rmse = c(res$rmse, rmse(tmp$AGE, tmp$single))
res$mae = c(res$mae, mae(tmp$AGE, tmp$single))
res$r = c(res$r, cor(tmp$AGE, tmp$single))
res$clock_used = c(res$clock_used, "single_cpg")
# # Single CpG no promoters
# my_coefs = feats_single_nopro %>%
#   dplyr::rename("Coef" = "s0")
# rownames(my_coefs) = my_coefs$IlmnID
# meta_whi$single_nopro = my_preds = predict_from_coefs(my_coefs, data_whi, untransform_age = 20)
# tmp = meta_whi %>%
#   drop_na(single_nopro)
# res$rmse = c(res$rmse, rmse(tmp$AGE, tmp$single_nopro))
# res$mae = c(res$mae, mae(tmp$AGE, tmp$single_nopro))
# res$r = c(res$r, cor(tmp$AGE, tmp$single_nopro))
# res$clock_used = c(res$clock_used, "single_cpg_nopro")
# Aggregated
my_coefs = feats_aggreg %>%
  dplyr::rename("Coef" = "s0")
rownames(my_coefs) = my_coefs$reID
meta_whi$aggreg = predict_from_coefs(my_coefs, data_re2_whi, untransform_age = 20)
tmp = meta_whi %>%
  drop_na(aggreg)
res$rmse = c(res$rmse, rmse(tmp$AGE, tmp$aggreg))
res$mae = c(res$mae, mae(tmp$AGE, tmp$aggreg))
res$r = c(res$r, cor(tmp$AGE, tmp$aggreg))
res$clock_used = c(res$clock_used, "aggreg")

# plot(tmp$AGE, tmp$horvath_sb)
# plot(tmp$AGE, tmp$single)
# plot(tmp$AGE, tmp$horvath_pheno)
# plot(tmp$AGE, tmp$aggreg)
all.equal(rownames(meta_whi), rownames(data_re2_whi))

res = as.data.frame(res)
res$clock_used = factor(res$clock_used, levels=res$clock_used)
p1 = ggplot(res, aes(x=clock_used, y=rmse, fill=clock_used))+
  geom_bar(stat="identity")+
  scale_y_continuous(expand=expansion(mult=c(0, 0.1)))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  ggtitle("RMSE")
p2 = ggplot(res, aes(x=clock_used, y=mae, fill=clock_used))+
  geom_bar(stat="identity")+
  scale_y_continuous(expand=expansion(mult=c(0, 0.1)))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  ggtitle("MAE")
p3 = ggplot(res, aes(x=clock_used, y=r, fill=clock_used))+
  geom_bar(stat="identity")+
  scale_y_continuous(expand=expansion(mult=c(0, 0.1)))+
  coord_cartesian(ylim = c(0, 1))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  ggtitle("r")
p1+p2+p3+plot_layout(guides="collect")+plot_annotation(title="Test on WHI BA23")
ggsave(paste0(paths$out, "/comparison_on_WHI.png"), width=1*w, height=0.4*h, units="mm", dpi = 300)

##### Summary #####

clocks_oi = c("horvath_pt", "horvath_sb", "horvath_pheno", "single", "aggreg") #,"single_nopro"
meta_whi[, paste0(clocks_oi, "Accel")] = meta_whi[, clocks_oi] - meta_whi$AGE
meta_whi = meta_whi %>%
  mutate(statusMortality=ifelse(!is.na(DEATHDY) & DEATHDY <= ENDFOLLOWDY, 2, 1)) %>%
  mutate(statusCancer=ifelse(!is.na(ANYCANCERDY) & ANYCANCERDY <= ENDFOLLOWDY, 2, 1)) %>%
  mutate(statusCHD=ifelse(!is.na(CHDDY) & CHDDY <= ENDFOLLOWDY, 2, 1))

res_summary = data.frame()
for (cl in clocks_oi) {
  accel = paste0(cl, "Accel")
  # Mortality
  index = paste0(cl, "Mortality")
  frml = sprintf("Surv(ENDFOLLOWDY, statusMortality) ~ %s+AGE", accel)
  test = summary(coxph(as.formula(frml), data = meta_whi))
  res_summary[index, "clock"] = cl
  res_summary[index, "predicted"] = "Mortality"
  res_summary[index, "coef"] = test$coefficients[1, 1]
  res_summary[index, "se"] = test$coefficients[1, 3]
  res_summary[index, "pval"] = test$coefficients[1, 5]
  # Cancer
  index = paste0(cl, "Cancer")
  frml = sprintf("Surv(ENDFOLLOWDY, statusCancer) ~ %s+AGE", accel)
  test = summary(coxph(as.formula(frml), data = meta_whi))
  res_summary[index, "clock"] = cl
  res_summary[index, "predicted"] = "Cancer"
  res_summary[index, "coef"] = test$coefficients[1, 1]
  res_summary[index, "se"] = test$coefficients[1, 3]
  res_summary[index, "pval"] = test$coefficients[1, 5]
  # CHD
  index = paste0(cl, "CHD")
  frml = sprintf("Surv(ENDFOLLOWDY, statusCHD) ~ %s+AGE", accel)
  test = summary(coxph(as.formula(frml), data = meta_whi))
  res_summary[index, "clock"] = cl
  res_summary[index, "predicted"] = "CHD"
  res_summary[index, "coef"] = test$coefficients[1, 1]
  res_summary[index, "se"] = test$coefficients[1, 3]
  res_summary[index, "pval"] = test$coefficients[1, 5]
}

fig$M4$clock_risk_association = res_summary %>%
  mutate(clock = factor(clock, levels=clocks_oi[c(3,1,2,4,5,6)])) %>%
  mutate(sig = stars.pval(pval)) %>%
  ggplot(., aes(clock, coef, ymin=coef-se, ymax=coef+se, fill=clock, label=sig))+
  geom_bar(stat="identity")+
  geom_errorbar(width=0.2)+
  geom_text(aes(y=coef+se+0.005))+
  facet_wrap(~predicted)
ggsave(plot=fig$M4$clock_risk_association, paste0(paths$out, "/clock_risk_association.png"), width=1*w, height=0.4*h, units="mm", dpi = 300)

#### L1HS FOR BIO AGE ####

all.equal(rownames(meta_whi), rownames(data_re2_whi))
# Get average methylation of young L1s into meta_whi
which_fam = c("L1HS", "L1PA2", "L1PA3", "L1PA4")
which_fam = c("L1HS", "L1PA2", "L1PA3", "L1PA4", "L1MEi", "L1PA11", "L1MA4A", "L1M7")
which_id = re_info_whi1 %>%
  dplyr::filter(fam %in% which_fam) %>%
  mutate(fam = factor(fam, levels=which_fam)) %>%
  arrange(fam)
which_id = rownames(which_id)
meta_whi[, which_fam] = data_re1_whi[, which_id]
# Add indicator vars
meta_whi = meta_whi %>%
  mutate(statusMortality=ifelse(!is.na(DEATHDY) & DEATHDY <= ENDFOLLOWDY, 2, 1)) %>%
  mutate(statusCancer=ifelse(!is.na(ANYCANCERDY) & ANYCANCERDY <= ENDFOLLOWDY, 2, 1)) %>%
  mutate(statusCHD=ifelse(!is.na(CHDDY) & CHDDY <= ENDFOLLOWDY, 2, 1)) %>%
  mutate(cancer_in2y = statusCancer == 2 & ANYCANCERDY < 2*365) %>%
  mutate(cancer_in3y = statusCancer == 2 & ANYCANCERDY < 3*365) %>%
  mutate(cancer_in4y = statusCancer == 2 & ANYCANCERDY < 4*365) %>%
  mutate(cancer_in5y = statusCancer == 2 & ANYCANCERDY < 5*365)

# Does methylation at these elements decrease with age in these datasets
ggplot(meta_whi, aes(AGE, L1HS))+
  geom_point()+
  geom_smooth()

table(meta_whi$cancer_in2y)
table(meta_whi$cancer_in3y)
table(meta_whi$cancer_in4y)
table(meta_whi$cancer_in5y)

# test_and_box = function(data, xvar, yvar, age_correction=F) {
#   frml = sprintf("%s~AGE+%s", yvar, xvar)
#   test=summary(lm(as.formula(frml), data=data, na.action = "na.exclude"))
#   # print(test)
#   data$age_subtr = data[[yvar]] - test$coefficients[2,1] * data$AGE 
#   pval = test$coefficients[3, 4]
#   pval = ifelse(pval < 0.001, "p < 0.001", paste0("p = ", round(pval, 3)))
#   if (age_correction) {
#     p = ggplot(data=data)+
#       geom_boxplot(aes(x=.data[[xvar]], y=age_subtr, fill=.data[[xvar]]), outlier.shape = NA)+
#       annotate("text", x = 1.5, y=quantile(data[[yvar]], 0.99, na.rm=T), label=pval)+
#       labs(y=yvar)
#   } else {
#     p = ggplot(data=data)+
#       geom_boxplot(aes(x=.data[[xvar]], y=.data[[yvar]], fill=.data[[xvar]]), outlier.shape = NA)+
#       annotate("text", x = 1.5, y=quantile(data[[yvar]], 0.99, na.rm=T), label=pval)
#   }
#   return(list(p=p, test=test))
# }

# Cancer and L1 methylation
# test1 = test_and_box(meta_whi, "cancer_in3y", "L1HS", age_correction = T)$test
# test2 = test_and_box(meta_whi, "cancer_in3y", "L1M7", age_correction = T)$test
# print(test1$coefficients)
# print(test2)
# 
# ps = list()
# ps[[1]] = test_and_box(meta_whi, "cancer_in3y", "L1HS", age_correction = T)$p+coord_cartesian(ylim=c(0.79, 0.9))+ggtitle("L1HS")
# ps[[2]] = test_and_box(meta_whi, "cancer_in3y", "L1PA2", age_correction = T)$p+coord_cartesian(ylim=c(0.78, 0.88))+ggtitle("L1PA2")
# ps[[3]] = test_and_box(meta_whi, "cancer_in3y", "L1PA3", age_correction = T)$p+coord_cartesian(ylim=c(0.77, 0.87))+ggtitle("L1PA3")
# ps[[4]] = test_and_box(meta_whi, "cancer_in3y", "L1PA4", age_correction = T)$p+coord_cartesian(ylim=c(0.73, 0.86))+ggtitle("L1PA4")
# # ps[[5]] = test_and_box(meta_whi, "cancer_in3y", "L1MEi", age_correction = T)$p+coord_cartesian(ylim=c(0.64, 0.76))+ggtitle("L1MEi")
# # ps[[6]] = test_and_box(meta_whi, "cancer_in3y", "L1PA11", age_correction = T)$p+coord_cartesian(ylim=c(0.72, 0.85))+ggtitle("L1PA11")
# # ps[[7]] = test_and_box(meta_whi, "cancer_in3y", "L1MA4A", age_correction = T)$p+coord_cartesian(ylim=c(0.72, 0.84))+ggtitle("L1MA4A")
# # ps[[8]] = test_and_box(meta_whi, "cancer_in3y", "L1M7", age_correction = T)$p+coord_cartesian(ylim=c(0.65, 0.78))+ggtitle("L1M7")
# figM3$B = wrap_plots(ps)+plot_layout(nrow=1, guides="collect")
# ggsave(plot = figM3$B, paste0(paths$out, "/batch_mean_methylation.png"), width=1*w, height=0.3*h, units="mm", dpi = 300)

# Association with cancer and age at the same time
res = data.frame()
for (i in 1:length(which_fam)){
  frml = sprintf("%s~AGE+%s", which_fam[i], "cancer_in3y")
  test=summary(lm(as.formula(frml), data=meta_whi, na.action = "na.exclude"))
  res[i, "fam"] = which_fam[i]
  res[i, "coef_age"] = test$coefficients[2, 1]
  res[i, "pval_age"] = test$coefficients[2, 4]
  res[i, "coef_cancer"] = test$coefficients[3, 1]
  res[i, "pval_cancer"] = test$coefficients[3, 4]
  # # Subtract age effect
  # frml = sprintf("%s~AGE", which_fam[i])
  # fit = lm(as.formula(frml), data=meta_whi, na.action = "na.exclude")
  # meta_whi[paste0(which_fam[i], "_residual")] = meta_whi[which_fam[i]] - coef(fit)[2] * meta_whi$AGE
}
res = res %>%
  pivot_longer(cols=-c(fam), names_sep = "_", names_to=c("type", "var")) %>%
  pivot_wider(names_from = "type", values_from = "value") %>%
  mutate(fam = factor(fam, levels=which_fam)) %>%
  mutate(sig = ifelse(pval < 0.05, "*", ""))

maxp = max(-log10(res$pval))

p1 = res %>%
  dplyr::filter(var == "cancer") %>%
  ggplot(., aes(fam, var, size = -log10(pval), fill=coef, label=sig))+
  geom_point(pch=21)+
  geom_text(color="white")+
  theme(legend.position="top",
        axis.title = element_blank(),
        axis.text.x = element_blank())+
  scale_size(limits=c(0, maxp), breaks = c(0.5, 1, 1.5, 2))
p2 = res %>%
  dplyr::filter(var == "age") %>%
  ggplot(., aes(fam, var, size = -log10(pval), fill=coef, label=sig))+
  geom_point(pch=21)+
  geom_text(color="white")+
  theme(legend.position="bottom",
        axis.title.y = element_blank())+
  scale_size(limits=c(0, maxp), breaks = c(0.5, 1, 1.5, 2))
fig$M4$l1_cancer_assoc = p1
fig$M4$l1_age_assoc = p2
p1/p2

#### WHI CANCER DETAIL ####

# Find cancer types within 3y
meta_whi_cancer3y = meta_whi %>%
  dplyr::filter(cancer_in3y)
colnames(meta_whi_cancer3y)[grepl("CANCER", colnames(meta_whi_cancer3y))]
# Columns containing cancer data (from WHI data dictionary)
cancer_cols = c("SINUS", "ADRENAL", "ANAL", "APPENDIX", "BASETONGUE", "BILIARY", 
                "BLADDER", "BONELIMB", "BONENON", "BRAIN", "CNSCA", "CERVICAL", 
                "COLON", "COLORECTAL", "CONNECTIVE", "ENDOCRINE", "ENDMTRL", 
                "ESOPHAGUS", "EYE", "MOUTHFLOOR", "GALLBLDR", "GENITAL", "GUM", 
                "HEART", "HYPOPHAR", "KIDNEY", "LARYNX", "LEUKEMIA", "LIVER", 
                "LUNG", "LYMPH", "HODGKINS", "LYMPHOMA", "MELANOMA", "MENINGES", 
                "MMYELOMA", "MYCOSISF","NASAL", "NASOPHAR", "ORALUNSP", "OROPHARYNX", 
                "OTHERDIGEST", "OTHERLIP", "OVARY", "PALATE", "PANCREAS", 
                "PAROTID", "PERIPHERAL", "PERITONEUM", "PYRIFORM", "RECTOSIG", 
                "RECTUM", "RENALPELV", "RESP", "SALIVARY", "SMINTEST", "STOMACH", 
                "THYMUS", "THYROID", "TONGUE", "TONSIL", "TRACHEA", "URETER", 
                "URINARY", "UTERINE", "VAGINA", "VULVA", "CANCOTHER")

# Which of these are actually in the data?
tmp = colSums(meta_whi_cancer3y[, cancer_cols])
tmp[tmp > 0]
cancer_cols_found = names(tmp[tmp > 0])

# Are any of these redundant?
View(meta_whi_cancer3y[, cancer_cols_found])
rowSums(meta_whi_cancer3y[, cancer_cols_found])
cor(meta_whi_cancer3y[, cancer_cols_found])
# Colorectal encompasses colon, rectum and rectosig
cancer_cols_found = setdiff(cancer_cols_found, c("COLON", "RECTUM", "RECTOSIG"))
View(meta_whi_cancer3y[, cancer_cols_found])
rowSums(meta_whi_cancer3y[, cancer_cols_found])
cor(meta_whi_cancer3y[, cancer_cols_found])
# Combine bladder and ureter into urinary (only case is both anyways)
meta_whi_cancer3y$URINARY = meta_whi_cancer3y$BLADDER
cancer_cols_found = c(cancer_cols_found, "URINARY")
cancer_cols_found = setdiff(cancer_cols_found, c("BLADDER", "URETER"))
rowSums(meta_whi_cancer3y[, cancer_cols_found])
cor(meta_whi_cancer3y[, cancer_cols_found])

# Breast cancer data is found in another file
unk_patients = meta_whi_cancer3y[rowSums(meta_whi_cancer3y[, cancer_cols_found]) == 0, "SubjectID"]
bc_patients = fread("data/whi/meta/outc_bc_inv.dat")
all(unk_patients %in% as.character(bc_patients$ID))
meta_whi_cancer3y$BREAST = as.numeric(meta_whi_cancer3y$SubjectID %in% bc_patients$ID)
sum(meta_whi_cancer3y$BREAST)
cancer_cols_found = c(cancer_cols_found, "BREAST")
View(meta_whi_cancer3y[, cancer_cols_found])
rowSums(meta_whi_cancer3y[, cancer_cols_found])

nocancer_mean = meta_whi %>%
  dplyr::filter(!cancer_in3y) %>%
  dplyr::select(L1HS, L1PA2, L1PA3, L1PA4) %>%
  colMeans(na.rm = T)
cancer_mean = meta_whi %>%
  dplyr::filter(cancer_in3y) %>%
  dplyr::select(L1HS, L1PA2, L1PA3, L1PA4) %>%
  colMeans(na.rm = T)
nocancer_mean
cancer_mean

p1 = meta_whi_cancer3y %>%
  pivot_longer(cols = cancer_cols_found, names_to = "CancerType") %>%
  dplyr::filter(value == 1) %>%
  ggplot(., aes(CancerType, L1HS))+
  geom_point()+
  stat_summary(
    geom = "point",
    fun.y = "mean",
    col = "black",
    size = 3,
    shape = 24,
    fill = "red")+
  geom_hline(yintercept = nocancer_mean["L1HS"], linetype="dashed", color="blue")+
  geom_hline(yintercept = cancer_mean["L1HS"], linetype="dashed", color="red")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  ggtitle("L1HS")+
  labs(y="Average L1HS CpG methylation")

p2 = meta_whi_cancer3y %>%
  pivot_longer(cols = cancer_cols_found, names_to = "CancerType") %>%
  dplyr::filter(value == 1) %>%
  ggplot(., aes(CancerType, L1PA2))+
  geom_point()+
  stat_summary(
    geom = "point",
    fun.y = "mean",
    col = "black",
    size = 3,
    shape = 24,
    fill = "red")+
  geom_hline(yintercept = nocancer_mean["L1PA2"], linetype="dashed", color="blue")+
  geom_hline(yintercept = cancer_mean["L1PA2"], linetype="dashed", color="red")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  ggtitle("L1PA2")+
  labs(y="Average L1PA2 CpG methylation")

p3 = meta_whi_cancer3y %>%
  pivot_longer(cols = cancer_cols_found, names_to = "CancerType") %>%
  dplyr::filter(value == 1) %>%
  ggplot(., aes(CancerType, L1PA3))+
  geom_point()+
  stat_summary(
    geom = "point",
    fun.y = "mean",
    col = "black",
    size = 3,
    shape = 24,
    fill = "red")+
  geom_hline(yintercept = nocancer_mean["L1PA3"], linetype="dashed", color="blue")+
  geom_hline(yintercept = cancer_mean["L1PA3"], linetype="dashed", color="red")+
  theme(axis.text.x=element_text(angle=90, hjust=1),
        axis.title.x = element_blank())+
  ggtitle("L1PA3")+
  labs(y="Average L1PA3 CpG methylation")

p4 = meta_whi_cancer3y %>%
  pivot_longer(cols = cancer_cols_found, names_to = "CancerType") %>%
  dplyr::filter(value == 1) %>%
  ggplot(., aes(CancerType, L1PA4))+
  geom_point()+
  stat_summary(
    geom = "point",
    fun.y = "mean",
    col = "black",
    size = 3,
    shape = 24,
    fill = "red")+
  geom_hline(yintercept = nocancer_mean["L1PA4"], linetype="dashed", color="blue")+
  geom_hline(yintercept = cancer_mean["L1PA4"], linetype="dashed", color="red")+
  theme(axis.text.x=element_text(angle=90, hjust=1),
        axis.title.x = element_blank())+
  ggtitle("L1PA4")+
  labs(y="Average L1PA4 CpG methylation")

fig$S6$canctype_l1hs = p1
fig$S6$canctype_l1pa2 = p2
fig$S6$canctype_l1pa3 = p3
fig$S6$canctype_l1pa4 = p4
fig$S6$canctype_assembled = p1 + p2 + p3 + p4 + plot_layout(nrow=2)
fig$S6$canctype_assembled

#### DISEASE PREDICTOR ####

##### Prep #####

# # Labels
# all.equal(as.character(meta_whi$SampleID), rownames(data_whi))
# y = data.frame(
#   time = meta_whi$ANYCANCERDY,
#   status =  meta_whi$statusCancer - 1)
# # table(y$status, is.na(y$time))
# y[y$status == 0, "time"] = meta_whi[y$status == 0, "ENDFOLLOWDY"]

probe_info_yl1 = probe_info %>%
  dplyr::filter(superf == "LINE/L1") %>%
  dplyr::filter(fam %in% c("L1HS", "L1PA2", "L1PA3", "L1PA4")) %>%
  dplyr::filter(IlmnID %in% colnames(data_whi))
cpgs_yl1 = probe_info_yl1$IlmnID

probe_info_the1 = probe_info %>%
  dplyr::filter(superf == "LTR/ERVL-MaLR") %>%
  dplyr::filter(fam %in% c("THE1A", "THE1C")) %>%
  dplyr::filter(IlmnID %in% colnames(data_whi))
cpgs_the1 = probe_info_the1$IlmnID

X_yl1 = data_whi[, cpgs_yl1]
sum(is.na(X_yl1)) / (nrow(X_yl1) * ncol(X_yl1))
X_yl1 = as.matrix(X_yl1)

X_yl1_the1 = data_whi[, c(cpgs_yl1, cpgs_the1)]
sum(is.na(X_yl1_the1)) / (nrow(X_yl1_the1) * ncol(X_yl1_the1))
X_yl1_the1 = as.matrix(X_yl1_the1)

meta_whi = meta_whi %>%
  mutate(mortality_in3y = statusMortality == 2 & DEATHDY < 3*365) %>%
  mutate(chd_in3y = statusCHD == 2 & CHDDY < 3*365)

##### Train #####

disease_predictors = list()
load(file=paste0(paths$out, "/disease_predictors.Rdata"))

set.seed(1337)
disease_predictors[["cancer"]] = cv.glmnet(X_yl1, meta_whi$cancer_in3y, family = "binomial", type.measure = "auc", alpha=0.9)
set.seed(1337)
disease_predictors[["mortality"]] = cv.glmnet(X_yl1, meta_whi$mortality_in3y, family = "binomial", type.measure = "auc", alpha=0.9)
set.seed(1337)
disease_predictors[["chd"]] = cv.glmnet(X_yl1, meta_whi$chd_in3y, family = "binomial", type.measure = "auc", alpha=0.9)
set.seed(1337)
disease_predictors[["cancer_w_the1"]] = cv.glmnet(X_yl1_the1, meta_whi$cancer_in3y, family = "binomial", type.measure = "auc", alpha=0.9)
set.seed(1337)
disease_predictors[["mortality_w_the1"]] = cv.glmnet(X_yl1_the1, meta_whi$mortality_in3y, family = "binomial", type.measure = "auc", alpha=0.9)
set.seed(1337)
disease_predictors[["chd_w_the1"]] = cv.glmnet(X_yl1_the1, meta_whi$chd_in3y, family = "binomial", type.measure = "auc", alpha=0.9)

disease_predictors[["cancer"]]
disease_predictors[["cancer_w_the1"]]
disease_predictors[["mortality"]] 
disease_predictors[["mortality_w_the1"]]
disease_predictors[["chd"]]
disease_predictors[["chd_w_the1"]]

res = data.frame(
  predicted = rep(c("cancer", "mortality", "chd"), each=2),
  lambda = rep(c("1sd", "min"), n=3))
for(i in 1:3) {
  imin = disease_predictors[1:3][[i]]$index[1]
  i1sd = disease_predictors[1:3][[i]]$index[2]
  res[2*(i-1)+1, "mean"] = disease_predictors[1:3][[i]]$cvm[i1sd]
  res[2*(i-1)+1, "lower"] = disease_predictors[1:3][[i]]$cvlo[i1sd]
  res[2*(i-1)+1, "upper"] = disease_predictors[1:3][[i]]$cvup[i1sd]
  res[2*(i-1)+1, "n"] = disease_predictors[1:3][[i]]$nzero[i1sd]
  res[2*(i-1)+2, "mean"] = disease_predictors[1:3][[i]]$cvm[imin]
  res[2*(i-1)+2, "lower"] = disease_predictors[1:3][[i]]$cvlo[imin]
  res[2*(i-1)+2, "upper"] = disease_predictors[1:3][[i]]$cvup[imin]
  res[2*(i-1)+2, "n"] = disease_predictors[1:3][[i]]$nzero[imin]
}

res_w_the1 = data.frame(
  predicted = rep(c("cancer_w_the1", "mortality_w_the1", "chd_w_the1"), each=2),
  lambda = rep(c("1sd", "min"), n=3))
for(i in 1:3) {
  imin = disease_predictors[4:6][[i]]$index[1]
  i1sd = disease_predictors[4:6][[i]]$index[2]
  res_w_the1[2*(i-1)+1, "mean"] = disease_predictors[4:6][[i]]$cvm[i1sd]
  res_w_the1[2*(i-1)+1, "lower"] = disease_predictors[4:6][[i]]$cvlo[i1sd]
  res_w_the1[2*(i-1)+1, "upper"] = disease_predictors[4:6][[i]]$cvup[i1sd]
  res_w_the1[2*(i-1)+1, "n"] = disease_predictors[4:6][[i]]$nzero[i1sd]
  res_w_the1[2*(i-1)+2, "mean"] = disease_predictors[4:6][[i]]$cvm[imin]
  res_w_the1[2*(i-1)+2, "lower"] = disease_predictors[4:6][[i]]$cvlo[imin]
  res_w_the1[2*(i-1)+2, "upper"] = disease_predictors[4:6][[i]]$cvup[imin]
  res_w_the1[2*(i-1)+2, "n"] = disease_predictors[4:6][[i]]$nzero[imin]
}


fig$M4$disease_predictors = ggplot(res, aes(lambda, mean, ymin=lower, ymax=upper, fill=lambda))+
  geom_bar(stat="identity")+
  geom_errorbar(width=0.2)+
  facet_wrap(~predicted)+
  coord_cartesian(ylim=c(0.5, 0.7))


fig$S6$disease_predictors_w_the1 = ggplot(res_w_the1, aes(lambda, mean, ymin=lower, ymax=upper, fill=lambda))+
  geom_bar(stat="identity")+
  geom_errorbar(width=0.2)+
  facet_wrap(~predicted)+
  coord_cartesian(ylim=c(0.5, 0.7))

save(disease_predictors, file=paste0(paths$out, "/disease_predictors.Rdata"))

table(meta_whi$cancer_in3y)
table(meta_whi$mortality_in3y)
table(meta_whi$chd_in3y)

#### SAVE ####

save(fig, file=paste0(paths$out, "/figures.Rdata"))

disease_features = list(
  cancer_min = data.frame(as.matrix(coef(disease_predictors$cancer, s="lambda.min"))) %>%
    dplyr::filter(s1 != 0) %>%
    rownames_to_column("IlmnID") %>%
    merge(., probe_info, by="IlmnID", all.x=T),
  cancer_1sd = data.frame(as.matrix(coef(disease_predictors$cancer, s="lambda.1se"))) %>%
    dplyr::filter(s1 != 0)%>%
    rownames_to_column("IlmnID") %>%
    merge(., probe_info, by="IlmnID", all.x=T),
  chd_min = data.frame(as.matrix(coef(disease_predictors$chd, s="lambda.min"))) %>%
    dplyr::filter(s1 != 0)%>%
    rownames_to_column("IlmnID") %>%
    merge(., probe_info, by="IlmnID", all.x=T),
  chd_1sd = data.frame(as.matrix(coef(disease_predictors$chd, s="lambda.1se"))) %>%
    dplyr::filter(s1 != 0)%>%
    rownames_to_column("IlmnID") %>%
    merge(., probe_info, by="IlmnID", all.x=T),
  mortality_min = data.frame(as.matrix(coef(disease_predictors$mortality, s="lambda.min"))) %>%
    dplyr::filter(s1 != 0)%>%
    rownames_to_column("IlmnID") %>%
    merge(., probe_info, by="IlmnID", all.x=T),
  mortality_1sd = data.frame(as.matrix(coef(disease_predictors$mortality, s="lambda.1se"))) %>%
    dplyr::filter(s1 != 0)%>%
    rownames_to_column("IlmnID") %>%
    merge(., probe_info, by="IlmnID", all.x=T),
  cancer_the1_min = data.frame(as.matrix(coef(disease_predictors$cancer_w_the1, s="lambda.min"))) %>%
    dplyr::filter(s1 != 0) %>%
    rownames_to_column("IlmnID") %>%
    merge(., probe_info, by="IlmnID", all.x=T),
  cancer_the1_1sd = data.frame(as.matrix(coef(disease_predictors$cancer_w_the1, s="lambda.1se"))) %>%
    dplyr::filter(s1 != 0)%>%
    rownames_to_column("IlmnID") %>%
    merge(., probe_info, by="IlmnID", all.x=T),
  chd_the1_min = data.frame(as.matrix(coef(disease_predictors$chd_w_the1, s="lambda.min"))) %>%
    dplyr::filter(s1 != 0)%>%
    rownames_to_column("IlmnID") %>%
    merge(., probe_info, by="IlmnID", all.x=T),
  chd_the1_1sd = data.frame(as.matrix(coef(disease_predictors$chd_w_the1, s="lambda.1se"))) %>%
    dplyr::filter(s1 != 0)%>%
    rownames_to_column("IlmnID") %>%
    merge(., probe_info, by="IlmnID", all.x=T),
  mortality_the1_min = data.frame(as.matrix(coef(disease_predictors$mortality_w_the1, s="lambda.min"))) %>%
    dplyr::filter(s1 != 0)%>%
    rownames_to_column("IlmnID") %>%
    merge(., probe_info, by="IlmnID", all.x=T),
  mortality_the1_1sd = data.frame(as.matrix(coef(disease_predictors$mortality_w_the1, s="lambda.1se"))) %>%
    dplyr::filter(s1 != 0)%>%
    rownames_to_column("IlmnID") %>%
    merge(., probe_info, by="IlmnID", all.x=T)
)

for (dpred in names(disease_features)){
  write.csv(disease_features[[dpred]], paste0(paths$out, "/supp/", dpred, ".csv"))
}
