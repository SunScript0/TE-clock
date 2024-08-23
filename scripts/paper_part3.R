library(tidyverse)
library(patchwork)
library(glmnet)

setwd("/scratch/fmorandi/internal/RE_clock/")

source("./scripts/utils.R")

paths = list()
paths$data = "./data"
paths$out = "./results/paper_part3"

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
main_res = c("LINE", "SINE", "LTR", "DNA", "tRNA", "Simple_repeat")

#### FIGURE STORAGE ####

figRRBS = list()

#### LOAD DATA ####

load(paste0(paths$data, "/processed/rrbs_data.Rdata"))
re_stats = read.table(paste0(paths$data, "/processed/re_cpg_counts.tsv"), header = T) %>%
  mutate(class = factor(class, levels=re_order))%>%
  mutate(superf = factor(superf, levels = unique(.$superf[base::order(.$class, .$superf)]))) %>%
  dplyr::arrange(class)

# Here chose if we'll work on high-coverage, or all data
met = data_rrbs$met_hc
cov = data_rrbs$cov_hc
# Are there zero-cov REs?
sum(rowSums(cov == 0) > 0) 
# Only 4 REs have some NA values, might as well remove
keep = rowSums(cov == 0) == 0
met = met[keep, ]
cov = cov[keep, ]
perc = met / cov
any(is.na(perc))

# Get approrpiate re_info
re_info = data_rrbs$re_info[rownames(perc), ]

# Note which re families are in the data
re_stats$in_data = re_stats$fam_id %in% re_info$fam_id
table(re_stats$class, re_stats$in_data)

# Simplify treatments for epi_memory
meta_rrbs = meta_rrbs %>%
  mutate(treatment=fct_recode(treatment,
                              "castration" = "old castration",
                              "castration control" = "old castration control",
                              "castration" = "young castrated",
                              "castration control" = "young control castrated"))

# Dual treatments
meta_rrbs = meta_rrbs %>%
  mutate(trt1 = gsub("-.*", "", treatment)) %>%
  mutate(trt2 = str_extract(treatment, "-(.*)", group=1))

# Take out cell culture experiments from petkovich
special_meta = subset(meta_rrbs, study=="petkovich" & tissue!="whole blood")
meta = meta_rrbs[!rownames(meta_rrbs) %in% rownames(special_meta), ]
met = met[, rownames(meta)]
cov = cov[, rownames(meta)]
perc = perc[, rownames(meta)]

#### OUTLIER REMOVAL ####

tmp = perc %>%
  drop_na() %>%
  t() %>%
  data.frame()
outlier = outlier_detection(tmp, meta, group="study", 
                            subgroups = c("tissue", "sex", "trt1"),
                            savepath = paste0(paths$out, "/outliers.pdf"))
meta$outlier = outlier

meta = subset(meta, !outlier)
met = met[, rownames(meta)]
cov = cov[, rownames(meta)]
perc = perc[, rownames(meta)]

##### EXPLORATION #####

table(meta[c("study", "tissue")])

# Global methylation stuff
hist(unlist(perc), breaks=50)
median(unlist(perc), na.rm=T)
all(colnames(perc) == rownames(meta))
meta$median_methylation = apply(perc, 2, median)
ggplot(meta, aes(x=study, y=median_methylation))+
  geom_boxplot()
ggsave(paste0(paths$out, "/median_study_methylation.png"), width=0.7*w, height=0.35*h, units="mm", dpi = 300)
ggplot(meta, aes(x=age, y=median_methylation))+
  geom_jitter(width = 1, size=0.5)
mean(meta$age)

quantile_norm = function(data, ref=NULL) {
  if (is.null(ref)) {
    data2 = apply(data, 2, sort)
    ref = rowMeans(data2)
  }
  for (i in 1:ncol(data)) {
    data[, i] = ref[rank(data[, i])]
  }
  return(data)
}
perc2 = quantile_norm(perc)
plot(unlist(perc[seq(1,900,10),seq(1,900,10)]), unlist(perc2[seq(1,900,10),seq(1,900,10)]), pch=".")
hist(unlist(perc), breaks=50)
hist(unlist(perc2), breaks=50)

# Copynumber, n CpG
summary(re_stats$n)
summary(re_stats$ncpg)

p1 = re_stats %>%
  dplyr::filter(class %in% main_res) %>%
  ggplot(., aes(x=class, y=log10(n)))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle=45, hjust=1))
p2 = re_stats %>%
  dplyr::filter(class %in% main_res) %>%
  ggplot(., aes(x=class, y=log10(ncpg)))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle=45, hjust=1))
p1 / p2
ggsave(paste0(paths$out, "/class_ncpg.png"), width=0.6*w, height=0.5*h, units="mm", dpi = 300)

re_stats %>%
  dplyr::filter(class %in% main_res) %>%
  ggplot(., aes(x=log10(n), y=log10(ncpg), color=class))+
  geom_point(size=0.1)+
  geom_abline()+
  stat_ellipse() +
  guides(color = guide_legend(override.aes = list(size = 2)))
ggsave(paste0(paths$out, "/class_n_vs_ncpg.png"), width=1*w, height=0.6*h, units="mm", dpi = 300)

res = data.frame(
  ncpg=seq(1, 10000, 10)
)
for (i in 1:nrow(res)) {
  res[i, "nfams"] = sum(re_stats$ncpg > res$ncpg[i])
}
ggplot(res, aes(log10(ncpg), log10(nfams)))+
  geom_line()+
  geom_vline(xintercept = c(2, 3))

# Relative coverage
re_info$meanCov = 100 * rowMeans(cov)
re_info$meanRelCov = rowMeans(100 * sweep(cov, 2, colSums(cov), "/"))
ggplot(re_info, aes(x=class, y=log10(meanRelCov)))+
  geom_boxplot()

pca = prcomp(t(perc[, grepl("petkovich", colnames(perc))]), center=T, scale.=T)
pca = merge(pca$x[, c("PC1", "PC2")], meta, by=0)
ggplot(pca, aes(x=PC1, y=PC2, color=treatment))+
  geom_point()

#### CLOCKS ####

##### Train #####

# Universal settings
all(colnames(perc) == rownames(meta))
meta$untreated = meta$trt1 == "untreated" | (meta$trt1 == "WT" & meta$trt2 == "standard")
alphas = c(0, 0.1, 0.2, 0.5, 0.8, 0.9, 1)

source("./scripts/utils.R")

clocks = list()
# # === Base, healthy ===
# healthy_ind = which(meta$untreated)
# clocks$base = ncv(t(perc)[healthy_ind, ], ageT(meta[healthy_ind, "age"], adult_age = 6), 
#                   cv_function = cv_glmnet, cv_args = list(alphas=alphas))
# predict_plot(clocks$base$preds, untransform_age = 6)+
#   ggtitle("Base clock, healthy individuals", 
#           subtitle = sprintf("Lambda = %.3f, Alpha = %.3f", 
#                              clocks$base$fit$lambda, 
#                              clocks$base$fit$alpha))
# ggsave(paste0(paths$out, "/clock_base_healthy.png"), width=0.6*w, height=0.3*h, units="mm", dpi = 300)

# # === Quantile ===
# clocks$quant = ncv(t(perc2)[healthy_ind, ], ageT(meta[healthy_ind, "age"], adult_age = 6),seed=clocks$base$seed,
#                   cv_function = cv_glmnet, cv_args = list(alphas=alphas))
# predict_plot(clocks$quant$preds, untransform_age = 6)+
#   ggtitle("Quantile norm clock, healthy", 
#           subtitle = sprintf("Lambda = %.3f, Alpha = %.3f", 
#                              clocks$quant$fit$lambda, 
#                              clocks$quant$fit$alpha))
# ggsave(paste0(paths$out, "/clock_quant_healthy.png"), width=0.6*w, height=0.3*h, units="mm", dpi = 300)

# # === High mean coverage ===
# all(rownames(re_info) == rownames(perc))
# th = -2.2
# ggplot(re_info, aes(log10(meanRelCov)))+
#   geom_histogram(bins=100)+
#   geom_vline(xintercept = th, color="red")
# re_info$highRelCov = log10(re_info$meanRelCov) > th
# ggplot(re_info, aes(log10(meanCov), fill=highRelCov))+
#   geom_histogram(bins=100)
# table(re_info[re_info$highRelCov, "class"]) / table(re_info$class)
# # Although one might think filtering by meanCov would remove almost all SimpleRepeats
# # I guess there is enough CpG island bias to keep them in
# inds = which(re_info$highRelCov)
# 
# clocks$filt_by_cov = ncv(t(perc)[healthy_ind, inds], ageT(meta[healthy_ind, "age"], adult_age = 6),seed=clocks$base$seed,
#                          cv_function = cv_glmnet, cv_args = list(alphas=alphas))
# predict_plot(clocks$filt_by_cov$preds, untransform_age = 6)+
#   ggtitle("Filtered by mean rel cov, healthy", 
#           subtitle = sprintf("Lambda = %.3f, Alpha = %.3f", 
#                              clocks$filt_by_cov$fit$lambda, 
#                              clocks$filt_by_cov$fit$alpha))
# ggsave(paste0(paths$out, "/clock_filt_by_cov.png"), width=0.6*w, height=0.3*h, units="mm", dpi = 300)

# # === Healthy, no thompson ===
# healthy_noth_ind = which(meta$untreated & meta$study != "thompson")
# 
# clocks$noth = ncv(t(perc)[healthy_noth_ind, ], ageT(meta[healthy_noth_ind, "age"], adult_age = 6),seed=clocks$base$seed,
#                          cv_function = cv_glmnet, cv_args = list(alphas=alphas))
# predict_plot(clocks$noth$preds, untransform_age = 6)+
#   ggtitle("No Thompson, healthy", 
#           subtitle = sprintf("Lambda = %.3f, Alpha = %.3f", 
#                              clocks$noth$fit$lambda, 
#                              clocks$noth$fit$alpha))
# 
# ggsave(paste0(paths$out, "/clock_noth.png"), width=0.6*w, height=0.3*h, units="mm", dpi = 300)

# === Healthy, no thompson, TE only ===
all(colnames(perc) == rownames(meta))
all(rownames(perc) == rownames(re_info))

healthy_noth_ind = which(meta$untreated & meta$study != "thompson")
te_ind = which(re_info$class %in% c("LINE", "SINE", "LTR", "DNA"))

table(meta[healthy_noth_ind, "study"], meta[healthy_noth_ind, "tissue"])

clocks$noth_te = ncv(t(perc)[healthy_noth_ind, te_ind], ageT(meta[healthy_noth_ind, "age"], adult_age = 6),
                  cv_function = cv_glmnet, cv_args = list(alphas=alphas))
predict_plot(clocks$noth_te$preds, meta=meta, color_var = "study", untransform_age = 6)+
  ggtitle("No Thompson, healthy", 
          subtitle = sprintf("Lambda = %.3f, Alpha = %.3f", 
                             clocks$noth_te$fit$lambda, 
                             clocks$noth_te$fit$alpha))

# Paper plot
tmp = clocks$noth_te$preds %>%
  mutate(age = ageTinv(age, 6) / 4.3) %>%
  mutate(preds = ageTinv(preds, 6) / 4.3) %>%
  merge(meta, ., by=0)
anno = sprintf("RMSE = %f\nMAE = %f\nr = %f", 
               rmse(tmp$age.y, tmp$preds),
               mae(tmp$age.y,tmp$preds),
               cor(tmp$age.y, tmp$preds))
xrange = range(tmp$age.y)
yrange = range(tmp$preds)
figRRBS$clock = ggplot(tmp, aes(x=age.y, y=preds, color=study))+
  geom_point(size=0.4)+
  labs(x="Age [months]", y="Predicted age [months]")+
  geom_abline(intercept=-2, slope=1, linetype="dashed")+
  geom_abline(intercept=0, slope=1)+
  geom_abline(intercept=+2, slope=1, linetype="dashed")+
  annotate("text", x = xrange[1], y=yrange[2], label=anno, size=3, vjust=1, hjust=0)+
  lims(x=xrange, y=yrange)

ggsave(paste0(paths$out, "/clock_noth_te.png"), width=0.6*w, height=0.3*h, units="mm", dpi = 300)

# # === Healthy, petkovich only, TE only ===
# all(colnames(perc) == rownames(meta))
# all(rownames(perc) == rownames(re_info))
# 
# healthy_petk_ind = which(meta$untreated & meta$study == "petkovich")
# te_ind = which(re_info$class %in% c("LINE", "SINE", "LTR", "DNA"))
# 
# clocks$petk_te = ncv(t(perc)[healthy_petk_ind, te_ind], ageT(meta[healthy_petk_ind, "age"], adult_age = 6),
#                      cv_function = cv_glmnet, cv_args = list(alphas=alphas))
# predict_plot(clocks$petk_te$preds, meta=meta, color_var = "study", untransform_age = 6)+
#   ggtitle("Petkovich, healthy", 
#           subtitle = sprintf("Lambda = %.3f, Alpha = %.3f", 
#                              clocks$petk_te$fit$lambda, 
#                              clocks$petk_te$fit$alpha))
# 
# ggsave(paste0(paths$out, "/clock_petk_te.png"), width=0.6*w, height=0.3*h, units="mm", dpi = 300)

##### Save #####

# saveRDS(clocks, file = paste0(paths$out, "/clocks.RDS"))
clocks = readRDS(paste0(paths$out, "/clocks.RDS"))

##### Selected features #####

a = data.frame(as.matrix(coef(clocks$noth_te$fit)))
features = data.frame(as.matrix(coef(clocks$noth_te$fit))) %>%
  dplyr::filter(s0 != 0) %>%
  merge(re_info, ., by=0) %>%
  arrange(class, superf, fam)
write.csv(features, "./results/sup_files/rrbs_coef.csv")
tmp3 = features %>%
  dplyr::filter(!is.na(class)) %>%
  group_by(class) %>%
  dplyr::summarize(nfeat=n()) %>%
  mutate(clock = "Combined CpG\nTE clock\n(Mouse RRBS)")
tmp1 = read.csv(paste0(paths$out, "/../paper_part2/class_composition_single.csv"), row.names = 1) %>%
  mutate(clock = "Individual CpG\nTE clock\n(Human array)")
tmp2 = read.csv(paste0(paths$out, "/../paper_part2/class_composition_aggreg.csv"), row.names = 1) %>%
  mutate(clock = "Combined CpG\nTE clock\n(Human array)")


p = do.call(rbind, list(tmp1, tmp2, tmp3)) %>%
  mutate(clock = factor(clock, levels = unique(clock))) %>%
  ggplot(., aes(clock, nfeat, fill=class))+
  geom_bar(stat = "identity", position="fill")+
  scale_y_continuous(expand = expansion(mult=0), label = scales::percent)+
  labs(y="% clock of features")
save(p, file=paste0(paths$out, "/class_composition.Rdata"))


##### Treatments #####

treated_noth_ind = which(!meta$untreated & meta$study != "thompson")
meta_treated = meta[treated_noth_ind, ]
data_treated = t(perc)[treated_noth_ind, te_ind]
meta_treated$pred = ageTinv(predict(clocks$noth_te$fit, newx=data_treated, s="lambda.min"), adult_age=6)

meta_treated = meta_treated %>%
  mutate(age = age / 4.3, pred = pred / 4.3)

table(meta_treated$study)

ggplot(meta_treated, aes(age, pred))+
  geom_point()+
  facet_wrap(~study, scales = "free")

ggplot(meta_treated, aes(treatment, pred-age))+
  geom_boxplot()+
  facet_wrap(~study, scales = "free")

meta_treated %>%
  dplyr::filter(study == "petkovich") %>%
  ggplot(., aes(age, pred, color=treatment))+
  geom_point()+
  geom_abline(intercept=-2, slope=1, linetype="dashed")+
  geom_abline(intercept=0, slope=1)+
  geom_abline(intercept=+2, slope=1, linetype="dashed")


figRRBS$treatments = meta_treated %>%
  dplyr::filter(study == "petkovich") %>%
  dplyr::filter(trt1 != "WT") %>%
  mutate(Experiment = ifelse(grepl("GHRKO", trt1), "GHRKO", "Snell Dwarf")) %>%
  mutate(Treatment = fct_recode(trt1, WT = "GHRKO WT", WT = "Snell Dwarf Control")) %>%
  mutate(Treatment = factor(Treatment, levels = c("WT", "GHRKO", "Snell Dwarf"))) %>%
  ggplot(., aes(Treatment, pred-age, fill=Treatment))+
  geom_boxplot()+
  facet_wrap(~Experiment, scale="free")

a = meta_treated %>%
  dplyr::filter(study == "petkovich") %>%
  dplyr::filter(trt1 != "WT") %>%
  mutate(Experiment = ifelse(grepl("GHRKO", trt1), "GHRKO", "Snell Dwarf")) %>%
  mutate(Treatment = fct_recode(trt1, WT = "GHRKO WT", WT = "Snell Dwarf Control")) %>%
  mutate(Treatment = factor(Treatment, levels = c("WT", "GHRKO", "Snell Dwarf"))) %>%
  mutate(accel = pred-age)

summary(lm(accel~Treatment, data=subset(a, Experiment == "GHRKO")))
wilcox.test(accel~Treatment, data=subset(a, Experiment == "GHRKO"))
summary(lm(accel~Treatment, data=subset(a, Experiment == "Snell Dwarf")))
wilcox.test(accel~Treatment, data=subset(a, Experiment == "Snell Dwarf"))

##### Save fig #####

figRRBS$clock
save(figRRBS, file=paste0(paths$out, "/figRRBS.Rdata"))

##### Compare #####

clocks$meer = read.table(paste0(paths$data, "/published_clocks/meer_coefs.txt"), header=T, sep="\t")
rownames(clocks$meer) = with(clocks$meer, c("(Intercept)", paste0(Chr[-1], ":", Coord[-1])))

th_single = readRDS(paste0(paths$data, "/rrbs_mouse/thompson/single.rds"))
rownames(th_single) = with(th_single, paste0(chr, ":", start))
th_single = th_single[, -c(1:3)]
th_met = th_single[, seq(1, ncol(th_single), 2)]
th_cov = th_single[, seq(1, ncol(th_single), 2)] + th_single[, seq(2, ncol(th_single), 2)]
th_perc = th_met / th_cov
colnames(th_perc) = str_replace(colnames(th_perc), "numCs_", "")
th_perc = data.frame(t(th_perc), check.names = F)

preds = predict_from_coefs(clocks$meer, th_perc)
which(rownames(clocks$meer) == "(Intercept)")
feats = rownames(clocks$meer)[-1]
preds = as.matrix(th_perc[, feats]) %*% clocks$meer[feats, "Coef"]
sort(setdiff(feats, colnames(th_perc)))
a = data.frame(sort(colnames(th_perc)))

a = merge(clocks$meer, th_single, by.x=c("Chr", "Coord"), by.y=c("chr", "start"))

