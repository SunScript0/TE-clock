library(tidyverse)
library(data.table)
library(limma)
library(wateRmelon)
library(Rfast)
library(MASS)
library(patchwork)
library(ggpubr)
library(scales)
library(mgcv)
library(ggdendro)
library(wesanderson)
library(gridExtra)
library(egg)
library(rstatix)
library(msa)
library(seqinr)
library(ape)
library(zoo)
library(edgeR)
library(BSgenome.Hsapiens.UCSC.hg19)

setwd("/scratch/fmorandi/internal/RE_clock/")

source("./scripts/utils.R")

paths = list()
paths$data = "./data"
paths$out = "./results/final"

#### PLOTTING SETTINGS ####

w = 210 - 25.4*2 # was 174
h = 297 - 25.4*2 # was 230
w_in = w/25.4
h_in = h/25.4

theme_set(theme_light(base_size = 8))
update_geom_defaults("text", list(size = 8*0.35))

re_order = c("LINE", "SINE", "LTR", "DNA", "DNA?", "Retroposon", "RC", 
             "tRNA", "srpRNA", "snRNA", "scRNA", "rRNA", 
             "Satellite", "Simple_repeat", "Low_complexity",
             "Unspecified", "Unknown")

#### FUNCTIONS ####

get_density = function(x, y, ...) {
  dens = MASS::kde2d(x, y, ...)
  ix = findInterval(x, dens$x)
  iy = findInterval(y, dens$y)
  ii = cbind(ix, iy)
  return(dens$z[ii])
}

#### LOAD DATA ####

data_infi = readRDS(paste0(paths$data, "/processed/infinium_data.rds"))
meta_infi = fread(paste0(paths$data, "/processed/infinium_meta.tsv"), data.table=F, sep="\t")
probe_info = readRDS(paste0(paths$data, "/processed/infinium_probes.rds"))
allRE_locs = read.table(paste0(paths$data, "/processed/re_locations.tsv"))

probe_info = probe_info %>%
  mutate(class = str_replace(class, "Satellite/centr","Satellite")) %>%
  mutate(class = factor(class, levels = re_order))
probe_info$superf = factor(probe_info$superf, levels = unique(probe_info$superf[base::order(probe_info$class, probe_info$superf)]))
rownames(probe_info) = probe_info$IlmnID

meta_whi = read.table(paste0(paths$data, "/processed/whi_meta.tsv"))

#### FIGURE STORAGE ####

load(paste0(paths$out, "/figures.Rdata"))
# fig = list()
# fig$M1 = list() # Dataset into, unadjusted LINE1, LTR
# fig$M2 = list() # Adjusted LINE1, LTR, motif analysis
# fig$M3 = list() # Clocks on chronological
# fig$M4 = list() # Clocks on disease
# fig$S1 = list() # CpG locations
# fig$S2 = list() # SINE, DNA unadjusted and adjusted
# fig$S3 = list() # ATAC and RNA
# fig$S4 = list() # Extra motif analysis
# fig$S5 = list() # Clock TE composition, overlap with existing clocks
# fig$S6 = list() # Cancer detail, extra disease predictors

#### DATASETS ####

tmp = rbind(
  meta_infi %>%
    dplyr::select(study, age),
  meta_whi %>%
    dplyr::rename(age = "AGE") %>%
    mutate(study = "WHI BA23") %>%
    dplyr::select(study, age)) %>%
  mutate(study = factor(study, levels = rev(c(unique(meta_infi$study), "WHI BA23"))))

# Compare age distributions
p1 = ggplot(meta_infi, aes(x="ALL", y=age))+
  geom_violin()+
  coord_flip()+
  theme(axis.title.y = element_blank())
p2 = ggplot(tmp, aes(x=study, y=age, fill=study))+
  geom_violin()+
  coord_flip()+
  theme(axis.title.y = element_blank(),
        legend.position="none")
fig$M1$dataset_ages = p1/p2 + plot_layout(heights = c(1,4))
ggsave(plot=fig$M1$dataset_ages, paste0(paths$out, "/dataset_ages.png"), width=1*w, height=0.5*h, units="mm", dpi = 300)

#### GENOMIC LOCATIONS ####

allRE_locs$group = "All REs"
probe_locs = probe_info %>%
  group_by(annotation) %>%
  summarize(n = n()) %>%
  mutate(group = "All probes")
probeRE_locs = probe_info %>%
  dplyr::filter(!is.na(superf)) %>%
  group_by(annotation) %>%
  summarize(n = n()) %>%
  mutate(group = "RE probes")
probeTE_locs = probe_info %>%
  dplyr::filter(class %in% c("LINE", "SINE", "LTR", "DNA", "Retroposon")) %>%
  group_by(annotation) %>%
  summarize(n = n()) %>%
  mutate(group = "TE probes")
fig$S1$genomic_locs = rbind(allRE_locs, probe_locs, probeRE_locs, probeTE_locs) %>%
  mutate(group = factor(group, levels=c("All REs", "All probes", "RE probes", "TE probes"))) %>%
  ggplot(., aes(x=group, y=n, fill=annotation))+
  geom_bar(stat="identity", position="fill") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  scale_y_continuous(expand=c(0,0))+
  labs(fill="Annotation")
ggsave(plot=fig$S1$genomic_locs, paste0(paths$out, "/genomic_locs.png"), width=1*w, height=0.4*h, units="mm", dpi = 300)

tmp1 = probe_info %>%
  dplyr::filter(!is.na(superf)) %>%
  group_by(class) %>%
  summarize(totn = n(), avg_cpg_dens = mean(CpG_density))
tmp2 = probe_info %>%
  dplyr::filter(!is.na(superf)) %>%
  group_by(class, annotation) %>%
  summarize(n = n())
tmp3 = probe_info %>%
  dplyr::filter(!is.na(superf)) %>%
  distinct(class, superf, fam) %>%
  group_by(class) %>%
  summarize(nfam = n())
probe_stats = merge(tmp1, tmp2, by="class")
probe_stats = merge(probe_stats, tmp3, by="class")

p1 = tmp1 %>%
  dplyr::filter(totn > 100) %>%
  ggplot(., aes(x=class, y=totn))+
  geom_bar(stat="identity", fill="#2255aa")+
  coord_flip()+
  scale_x_discrete(position = "top", limits=rev)+
  scale_y_reverse(expand=c(0,0),
                  breaks = breaks_pretty(n = 3))+
  labs(y="#CpGs")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
p2 = probe_stats %>%
  dplyr::filter(totn > 100) %>%
  ggplot(., aes(x=class, y=n, fill=annotation))+
  geom_bar(stat="identity", position="fill")+
  coord_flip()+
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position="bottom")+
  scale_y_continuous(expand=c(0,0))+
  scale_x_discrete(limits=rev)+
  labs(fill="Genomic location")+
  geom_text(data=subset(tmp1, totn > 100), aes(y=0, x=class, label=class, fill=NULL), hjust=0, color="white")
p3 = tmp1 %>%
  dplyr::filter(totn > 100) %>%
  ggplot(., aes(x=1, y=class, fill=avg_cpg_dens))+
  geom_tile(width=0.9, height=0.9)+
  scale_y_discrete(limits=rev)+
  scale_fill_viridis_c()+
  coord_fixed(ratio=1)+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  labs(fill="CpG density")

fig$S1$genomic_locs_class = p1+p2+p3+plot_layout(guides="keep", widths=c(4,16,1))
ggsave(plot = fig$S1$genomic_locs_class, paste0(paths$out, "/genomic_locs_class.png"), width=1*w, height=0.3*h, units="mm", dpi = 300)

rm(allRE_locs, probe_locs, probeRE_locs)

#### AGE COEFS ####

##### Combined #####

# all.equal(rownames(data_infi), meta_infi$ID)
# keep = meta_infi$study != "GSE87648" & meta_infi$health == "Healthy"
# design = model.matrix(~age+study+sex, data = meta_infi[keep, ])
# tmp = t(data_infi[keep, ])
# fit = lmFit(tmp, design)
# fit = eBayes(fit)
# coefs = coef(fit)
# res = topTable(fit, coef = "age", number = Inf, sort.by="none") # coef == logFC for age, in this case
# all.equal(rownames(coefs), rownames(res))
# all(coefs[,2] == res[1])
# res_infinium = data.frame(
#   age_coef = 100*res$logFC, # percentage rather than ratio, for nicer numbers
#   beta20 = 100*(coefs[,1] + coefs[,2]*20),
#   pval = res$P.Value,
#   padj = res$adj.P.Val)
# 
# save(res_infinium, file=paste0(paths$out, "/limma_results.Rdata"))

##### By study #####

# all.equal(rownames(data_infi),  meta_infi$ID)
# res_by_study = list()
# for (study in c("GSE40279", "GSE64495", "GSE147221", "GSE157131")) { # "GSE87648" excluded
#   keep = meta_infi$study == study & meta_infi$health == "Healthy"
#   meta_sset = meta_infi[keep, ]
#   data_sset = t(data_infi[meta_sset$ID, ])
#   design = model.matrix(~age+sex, data = meta_sset)
#   fit = lmFit(data_sset, design)
#   fit = eBayes(fit)
#   res = topTable(fit, coef = "age", number = Inf, sort.by="none")
#   coefs = coef(fit)
#   res_by_study[[study]] = data.frame(
#     IlmnID = rownames(res),
#     study = study,
#     age_coef = 100*res$logFC,
#     beta20 = 100*(coefs[,1] + coefs[,2]*20),
#     pval = res$P.Value,
#     padj = res$adj.P.Val,
#     se = sqrt(fit$s2.post) * fit$stdev.unscaled)
# }
# res_by_study = do.call(rbind, res_by_study)
# save(res_by_study, file=paste0(paths$out, "/limma_results_by_study.Rdata"))

##### Load #####

load(paste0(paths$out, "/limma_results.Rdata"))
# load(paste0("results/paper_part1/limma_results.Rdata")) # Identical as new
load(paste0(paths$out, "/limma_results_by_study.Rdata"))

probe_info = merge(probe_info, res_infinium, by.x="IlmnID", by.y=0)

##### Study reliability #####

all.equal(rep(rownames(res_infinium), times=4), res_by_study$IlmnID)
res_by_study2 = res_by_study %>%
  group_by(IlmnID) %>%
  summarize(age_coef = mean(age_coef), beta20 = mean(beta20))
all.equal(res_by_study2$IlmnID, rownames(res_infinium))
cor(res_by_study2$age_coef, res_infinium$age_coef)

res_by_study %>%
  dplyr::select(IlmnID, study, age_coef) %>%
  pivot_wider(names_from = "study", values_from="age_coef") %>%
  column_to_rownames("IlmnID") %>%
  cor() %>%
  as.data.frame() %>%
  rownames_to_column("Study1") %>%
  pivot_longer(cols=-c("Study1"), names_to = "Study2") %>%
  ggplot(., aes(Study1, Study2, fill=value, label=round(value, 2))) +
  geom_tile()+
  geom_text()+
  scale_fill_gradient2(low="blue", high="red")
ggsave(paste0(paths$out, "/meta_analysis_coef_cor.png"), width=0.6*w, height=0.3*h, units="mm", dpi = 300)

#### ADJUSTED COEFS ####

p1 = ggplot(probe_info, aes(x=beta20, y=age_coef))+
  geom_point(pch=".")+
  geom_smooth(se=F)+
  geom_hline(yintercept=0, color="red", linetype = "dashed")+
  labs(x="Methylation % at 20 yo", y="Methylation drift rate [%/year]")+
  ylim(-0.3, 0.3)
p2 = ggplot(probe_info, aes(x=CpG_density, y=age_coef))+
  geom_jitter(pch=".")+
  geom_smooth(se=F)+
  geom_hline(yintercept=0, color="red", linetype = "dashed")+
  labs(x="CpG density (100bp window)", y="Methylation drift rate [%/year]")+
  ylim(-0.3, 0.3)
fig$M2$adj_model1 = p1
fig$M2$adj_model2 = p2
ggsave(plot = p1+p2, paste0(paths$out, "/baseline_cpg_dens_vs_coef.png"), width=1*w, height=0.4*h, units="mm", dpi = 300)

probe_info %>%
  mutate(dens = get_density(CpG_density, beta20, n=100)) %>%
  ggplot(., aes(x=CpG_density, y=beta20, color=dens))+
  geom_jitter(pch=".")+
  scale_color_viridis_c()+
  labs(x="CpG density (100bp window)", y="Methylation % at 20 yo")+
  guides(color="none")
ggsave(paste0(paths$out, "/baseline_vs_cpg_dens.png"), width=0.6*w, height=0.3*h, units="mm", dpi = 300)

# Model to account for CpG density and baseline
model = gam(age_coef ~ s(beta20, bs = "cs") + s(CpG_density, bs = "cs") + s(beta20, bs = "cs", by=CpG_density), data=probe_info)
summary(model)

probe_info$residuals = residuals(model)

p1 = ggplot(probe_info, aes(x=beta20, y=residuals))+
  geom_point(pch=".")+
  geom_smooth(se=F)+
  geom_hline(yintercept=0, color="red", linetype = "dashed")+
  labs(x="Methylation % at 20 yo", y="Adj. methylation drift rate [%/year]")+
  ylim(-0.3, 0.3)
p2 = ggplot(probe_info, aes(x=CpG_density, y=residuals))+
  geom_jitter(pch=".")+
  geom_smooth(se=F)+
  geom_hline(yintercept=0, color="red", linetype = "dashed")+
  labs(x="CpG density (100bp window)", y="Adj. methylation drift rate [%/year]")+
  ylim(-0.3, 0.3)
p = p1+p2+plot_layout(nrow=1)
ggsave(plot = p, paste0(paths$out, "/baseline_cpg_dens_vs_coef_adj.png"), width=1*w, height=0.4*h, units="mm", dpi = 300)

#### AGE TRENDS ####

##### Global patterns #####

probe_info = probe_info %>%
  mutate(selfish = ifelse(class %in% c("LINE", "SINE", "LTR", "DNA", "Retroposon"), "Selfish", "Non-selfish"))

p1 = ggplot(probe_info, aes(x=beta20, fill=!is.na(superf)))+
  geom_density(alpha = 0.5)+
  labs(x="Methylation [% at 20 yo]", fill="In RE")
p2 = ggplot(probe_info, aes(x=age_coef, fill=!is.na(superf)))+
  geom_density(alpha=0.5, bw=0.004)+
  scale_y_sqrt()+
  xlim(-0.25, 0.25)+
  labs(x="Methylation drift rate [%/year]", fill="In RE")
fig$M1$all_cpg_baseline_and_coef = p1/p2+plot_layout(guides="collect")
ggsave(plot=fig$M1$all_cpg_baseline_and_coef, paste0(paths$out, "/all_cpg_baseline_and_coef.png"), width=0.6*w, height=0.5*h, units="mm", dpi = 300)

p1 = ggplot(probe_info, aes(x=beta20, fill=!is.na(superf)))+
  geom_density(alpha = 0.5)+
  labs(x="Methylation [% at 20 yo]", fill="In RE")
p2 = ggplot(probe_info, aes(x=residuals, fill=!is.na(superf)))+
  geom_density(alpha=0.5, bw=0.004)+
  scale_y_sqrt()+
  xlim(-0.25, 0.25)+
  labs(x="Methylation drift rate [%/year]", fill="In RE")
p1+p2+plot_layout(guides="collect")
ggsave(paste0(paths$out, "/all_cpg_baseline_and_coef_adj.png"), width=1*w, height=0.3*h, units="mm", dpi = 300)

##### By class #####

stats_class = probe_info %>%
  group_by(class) %>%
  summarize(med_age_coef = median(age_coef), med_res = median(residuals),
            mad_age_coef = mad(age_coef), mad_res = mad(residuals),
            mean_len = mean(re_len), n = n(), n_uniq = length(unique(rep_name))) %>%
  drop_na()
wilcox.test(age_coef ~ selfish, data=probe_info)
t.test(age_coef ~ selfish, data=probe_info)

fig$M1$class_age_coef = probe_info %>%
  dplyr::filter(!is.na(class)) %>%
  left_join(stats_class, by="class") %>%
  dplyr::filter(n > 100) %>%
  mutate(selfish = ifelse(class %in% c("LINE", "SINE", "LTR", "DNA", "Retroposon"), "Selfish", "Non-selfish")) %>%
  ggplot(., aes(x=class, y=age_coef, fill=selfish))+
  geom_boxplot(outlier.shape=NA)+
  geom_hline(yintercept=0, linetype="dashed", color="dark green") +
  scale_x_discrete(limits=rev)+
  scale_fill_manual(values=c("#1177dd", "#dd1111"))+
  coord_flip()+
  labs(y="Methylation drift rate [%/year]", fill="Repeat type")+
  ylim(-0.05, 0.05)+
  theme(axis.title.y=element_blank())+
  ggtitle("Methylation drift of RE classes", subtitle = "Classes with > 100 CpGs in array shown")
ggsave(plot = fig$M1$class_age_coef, paste0(paths$out, "/class_age_coef.png"), width=1*w, height=0.5*h, units="mm", dpi = 300)

probe_info %>%
  dplyr::filter(!is.na(class)) %>%
  left_join(stats_class, by="class") %>%
  dplyr::filter(n > 100) %>%
  mutate(selfish = ifelse(class %in% c("LINE", "SINE", "LTR", "DNA", "Retroposon"), "Selfish", "Non-selfish")) %>%
  ggplot(., aes(x=class, y=residuals, fill=selfish))+
  geom_boxplot(outlier.shape=NA)+
  geom_hline(yintercept=0, linetype="dashed", color="dark green") +
  scale_x_discrete(limits=rev)+
  scale_fill_manual(values=c("#1177dd", "#dd1111"))+
  coord_flip()+
  labs(y="Methylation drift rate [%/year]", fill="Repeat type")+
  ylim(-0.05, 0.05)+
  theme(axis.title.y=element_blank())+
  ggtitle("Methylation drift of RE classes", subtitle = "Classes with > 100 CpGs in array shown")
ggsave(paste0(paths$out, "/class_age_coef_adj.png"), width=1*w, height=0.5*h, units="mm", dpi = 300)

##### LINE/L1 #####

stats_l1 = probe_info %>%
  dplyr::filter(superf=="LINE/L1") %>%
  group_by(superf, fam) %>%
  summarize(med_age_coef = median(age_coef), med_res = median(residuals),
            mad_age_coef = mad(age_coef), mad_res = mad(residuals),
            mean_len = mean(re_len), n = n(), n_uniq = length(unique(rep_name))) %>%
  drop_na()

min_n=40

p1 = stats_l1 %>%
  dplyr::filter(n > min_n) %>%
  mutate(fam = fct_reorder(fam, mean_len)) %>%
  ggplot(., aes(x=fam, y=mean_len, group=1))+
  geom_line()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+
  labs(y="Mean seq\nlength")
p2 = probe_info %>%
  left_join(stats_l1, by="fam") %>%
  dplyr::filter(n > min_n) %>%
  mutate(rep = fct_reorder(fam, mean_len)) %>%
  ggplot(., aes(x=rep, y=age_coef, fill=after_stat(middle))) +
  geom_boxplot(outlier.shape=NA) +
  geom_hline(yintercept=0, color="red")+
  scale_fill_gradient2(low=muted("blue"), high=muted("red"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  labs(x="Subfamilies sorted by mean length", y="Methylation drift rate [%/year]",
       fill="Median methylation\ndrift rate [%/year]")+
  coord_cartesian(ylim=c(-0.11, 0.06))
fig$M1$l1_age_coef = p1/p2+plot_layout(heights = c(1,4))
ggsave(plot=fig$M1$l1_age_coef, paste0(paths$out, "/l1_age_coef.png"), width=1.5*w, height=0.4*h, units="mm", dpi = 300)

p2 = probe_info %>%
  dplyr::filter(superf=="LINE/L1") %>%
  left_join(stats_l1, by="fam") %>%
  dplyr::filter(n > min_n) %>%
  mutate(fam = fct_reorder(fam, mean_len)) %>%
  ggplot(., aes(x=fam, y=residuals, fill=after_stat(middle))) +
  geom_boxplot(outlier.shape=NA) +
  geom_hline(yintercept=0, color="red")+
  scale_fill_gradient2(low=muted("blue"), high=muted("red"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  labs(x="Subfamilies sorted by mean length", y="Adj. methylation drift rate [%/year]",
       fill="Median adj. methylation\ndrift rate [%/year]")+
  coord_cartesian(ylim=c(-0.07, 0.051))
fig$M2$l1_age_coef_adj = p1/p2+plot_layout(heights = c(1,4))
ggsave(plot = fig$M2$l1_age_coef_adj, paste0(paths$out, "/l1_age_coef_adj.png"), width=1.5*w, height=0.4*h, units="mm", dpi = 300)

##### SINE #####

stats_sine = probe_info %>%
  dplyr::filter(class=="SINE") %>%
  group_by(superf, fam) %>%
  summarize(med_age_coef = median(age_coef), med_res = median(residuals),
            mad_age_coef = mad(age_coef), mad_res = mad(residuals),
            mean_len = mean(re_len), n = n(), n_uniq = length(unique(rep_name))) %>%
  drop_na()

min_n=40

fig$S2$sine_age_coef = probe_info %>%
  dplyr::filter(class=="SINE") %>%
  left_join(stats_sine, by="fam") %>%
  dplyr::filter(n > min_n) %>%
  mutate(fam = fct_reorder(fam, mean_len)) %>%
  ggplot(., aes(x=fam, y=age_coef, fill=after_stat(middle))) +
  geom_boxplot(outlier.shape=NA)+
  geom_hline(yintercept=0, color="red")+
  scale_fill_gradient2(low=muted("blue"), high=muted("red"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  labs(x="Subfamilies sorted by mean length", y="Methylation drift rate [%/year]",
       fill="Median methylation\ndrift rate [%/year]")+
  coord_cartesian(ylim = c(-0.1, 0.06))
ggsave(plot = fig$S2$sine_age_coef, paste0(paths$out, "/sine_age_coef.png"), width=1*w, height=0.4*h, units="mm", dpi = 300)

fig$S2$sine_age_coef_adj = probe_info %>%
  dplyr::filter(class=="SINE") %>%
  left_join(stats_sine, by="fam") %>%
  dplyr::filter(n > min_n) %>%
  mutate(fam = fct_reorder(fam, mean_len)) %>%
  ggplot(., aes(x=fam, y=residuals, fill=after_stat(middle))) +
  geom_boxplot(outlier.shape=NA)+
  geom_hline(yintercept=0, color="red")+
  scale_fill_gradient2(low=muted("blue"), high=muted("red"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  labs(x="Subfamilies sorted by mean length", y="Adj. methylation drift rate [%/year]",
       fill="Median adj. methylation\ndrift rate [%/year]")+
  coord_cartesian(ylim = c(-0.06, 0.05))
ggsave(plot = fig$S2$sine_age_coef_adj, paste0(paths$out, "/sine_age_coef_adj.png"), width=1*w, height=0.4*h, units="mm", dpi = 300)

##### Retroposon/SVA #####

stats_sva = probe_info %>%
  dplyr::filter(class=="Retroposon") %>%
  group_by(superf, fam) %>%
  summarize(med_age_coef = median(age_coef), med_res = median(residuals),
            mad_age_coef = mad(age_coef), mad_res = mad(residuals),
            mean_len = mean(re_len), n = n(), n_uniq = length(unique(rep_name))) %>%
  drop_na()

min_n=40

fig$S2$sva_age_coef = probe_info %>%
  dplyr::filter(class=="Retroposon") %>%
  left_join(stats_sva, by="fam") %>%
  dplyr::filter(n > min_n) %>%
  mutate(fam = fct_reorder(fam, mean_len)) %>%
  ggplot(., aes(x=fam, y=age_coef, fill=after_stat(middle))) +
  geom_boxplot(outlier.shape=NA)+
  geom_hline(yintercept=0, color="red")+
  scale_fill_gradient2(low=muted("blue"), high=muted("red"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  labs(x="Subfamilies sorted by mean length", y="Methylation drift rate [%/year]",
       fill="Median methylation\ndrift rate [%/year]")+
  coord_cartesian(ylim = c(-0.05, 0.025))
ggsave(plot = fig$S2$sva_age_coef, paste0(paths$out, "/sva_age_coef.png"), width=1*w, height=0.4*h, units="mm", dpi = 300)

fig$S2$sva_age_coef_adj = probe_info %>%
  dplyr::filter(class=="Retroposon") %>%
  left_join(stats_sva, by="fam") %>%
  dplyr::filter(n > min_n) %>%
  mutate(fam = fct_reorder(fam, mean_len)) %>%
  ggplot(., aes(x=fam, y=residuals, fill=after_stat(middle))) +
  geom_boxplot(outlier.shape=NA)+
  geom_hline(yintercept=0, color="red")+
  scale_fill_gradient2(low=muted("blue"), high=muted("red"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  labs(x="Subfamilies sorted by mean length", y="Adj. methylation drift rate [%/year]",
       fill="Median adj. methylation\ndrift rate [%/year]")+
  coord_cartesian(ylim = c(-0.04, 0.025))
ggsave(plot = fig$S2$sva_age_coef_adj, paste0(paths$out, "/sva_age_coef_adj.png"), width=1*w, height=0.4*h, units="mm", dpi = 300)

##### LTR #####

stats_ltr = probe_info %>%
  dplyr::filter(class=="LTR") %>%
  group_by(superf, fam) %>%
  summarize(med_age_coef = median(age_coef), med_res = median(residuals),
            mad_age_coef = mad(age_coef), mad_res = mad(residuals),
            mean_len = mean(re_len), n = n(), n_uniq = length(unique(rep_name))) %>%
  drop_na()

min_n=40

fig$M1$ltr_age_coef = probe_info %>%
  ungroup() %>%
  dplyr::filter(class=="LTR") %>%
  left_join(stats_ltr, by=c("superf", "fam")) %>%
  dplyr::filter(n > min_n) %>%
  mutate(fam = fct_reorder(fam, mean_len)) %>%
  mutate(superf = str_remove(superf, "LTR/"))%>%
  ggplot(., aes(x=fam, y=age_coef, fill=after_stat(middle))) +
  geom_boxplot(outlier.shape=NA)+
  geom_hline(yintercept=0, color="red")+
  scale_fill_gradient2(low=muted("blue"), high=muted("red"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  labs(x="Subfamilies sorted by mean length", y="Methylation drift rate [%/year]",
       fill="Median methylation\ndrift rate [%/year]")+
  coord_cartesian(ylim = c(-0.12, 0.07))+
  facet_grid(~superf, scale="free_x", space = "free_x")
ggsave(plot=fig$M1$ltr_age_coef, paste0(paths$out, "/ltr_age_coef.png"), width=1.5*w, height=0.4*h, units="mm", dpi = 300)

fig$M2$ltr_age_coef_adj = probe_info %>%
  ungroup() %>%
  dplyr::filter(class=="LTR") %>%
  left_join(stats_ltr, by=c("superf", "fam")) %>%
  dplyr::filter(n > min_n) %>%
  mutate(fam = fct_reorder(fam, mean_len)) %>%
  mutate(superf = str_remove(superf, "LTR/"))%>%
  ggplot(., aes(x=fam, y=residuals, fill=after_stat(middle))) +
  geom_boxplot(outlier.shape=NA)+
  geom_hline(yintercept=0, color="red")+
  scale_fill_gradient2(low=muted("blue"), high=muted("red"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  labs(x="Subfamilies sorted by mean length", y="Adj. methylation drift rate [%/year]",
       fill="Median adj. methylation\ndrift rate [%/year]")+
  coord_cartesian(ylim = c(-0.09, 0.07))+
  facet_grid(~superf, scale="free_x", space = "free_x")
ggsave(plot =fig$M2$ltr_age_coef_adj, paste0(paths$out, "/ltr_age_coef_adj.png"), width=1.5*w, height=0.4*h, units="mm", dpi = 300)

##### DNA #####

stats_dna = probe_info %>%
  dplyr::filter(class=="DNA") %>%
  group_by(superf, fam) %>%
  summarize(med_age_coef = median(age_coef), med_res = median(residuals),
            mad_age_coef = mad(age_coef), mad_res = mad(residuals),
            mean_len = mean(re_len), n = n(), n_uniq = length(unique(rep_name))) %>%
  drop_na()

min_n=40

fig$S2$dna_age_coef = probe_info %>%
  dplyr::filter(class=="DNA") %>%
  left_join(stats_dna, by="fam") %>%
  dplyr::filter(n > min_n) %>%
  mutate(fam = fct_reorder(fam, mean_len)) %>%
  ggplot(., aes(x=fam, y=age_coef, fill=after_stat(middle))) +
  geom_boxplot(outlier.shape=NA)+
  geom_hline(yintercept=0, color="red")+
  scale_fill_gradient2(low=muted("blue"), high=muted("red"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  labs(x="Subfamilies sorted by mean length", y="Methylation drift rate [%/year]",
       fill="Median methylation\ndrift rate [%/year]")+
  coord_cartesian(ylim = c(-0.1, 0.06))
ggsave(plot = fig$S2$dna_age_coef, paste0(paths$out, "/dna_age_coef.png"), width=1*w, height=0.4*h, units="mm", dpi = 300)

fig$S2$dna_age_coef_adj = probe_info %>%
  dplyr::filter(class=="DNA") %>%
  left_join(stats_dna, by="fam") %>%
  dplyr::filter(n > min_n) %>%
  mutate(fam = fct_reorder(fam, mean_len)) %>%
  ggplot(., aes(x=fam, y=residuals, fill=after_stat(middle))) +
  geom_boxplot(outlier.shape=NA)+
  geom_hline(yintercept=0, color="red")+
  scale_fill_gradient2(low=muted("blue"), high=muted("red"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  labs(x="Subfamilies sorted by mean length", y="Adj. methylation drift rate [%/year]",
       fill="Median adj. methylation\ndrift rate [%/year]")+
  coord_cartesian(ylim = c(-0.065, 0.05))
ggsave(plot=fig$S2$dna_age_coef_adj, paste0(paths$out, "/dna_age_coef_adj.png"), width=1*w, height=0.4*h, units="mm", dpi = 300)

#### NON-LINEAR CHANGES ####

probes_l1 = probe_info %>%
  dplyr::filter(superf == "LINE/L1")
data_l1 = data_infi[,probes_l1$IlmnID]
data_l1_collapsed = collapse_to_clusters(data_l1, probes_l1$fam)

show_which = c("L1HS", "L1PA2", "L1PA3", "L1PA4", "L1MEi", "L1PA11", "L1MA4A", "L1M7")

all.equal(rownames(data_l1), meta_infi$ID)
tmp = meta_infi %>%
  cbind(., data_l1_collapsed[, show_which]) %>%
  dplyr::filter(study == "GSE40279") %>%
  # dplyr::filter(study != "GSE87648 ") %>%
  pivot_longer(cols = show_which, names_to="L1_fam") %>%
  mutate(L1_fam = factor(L1_fam, levels = show_which))
min(tmp$value)

fig$M4$non_linearity = tmp %>%
  dplyr::filter(value < 0.9) %>%
  ggplot(., aes(x=age, y=100*value)) +
  geom_point(size=0.2)+
  geom_smooth()+
  geom_smooth(data = subset(tmp, age < 65), method="lm", color="orange", linetype="dashed",fullrange=TRUE)+
  facet_wrap(~L1_fam, nrow=2, scales = "free_y")+
  labs(x="Age", y="Methylation [%]")
ggsave(plot = fig$M4$non_linearity, paste0(paths$out, "/non_linearity.png"), width=1*w, height=0.4*h, units="mm", dpi = 300)

#### MULTIOMIC INTERATION ####

##### RNA #####

load(paste0(paths$data, "/processed/rna_data.Rdata"))

# Fit EdgeR
dge_rna = DGEList(counts_re, samples = meta_rna)
design = model.matrix(~Age+Sex, data=meta_rna)
dge_rna = calcNormFactors(dge_rna)
dge_rna = estimateDisp(dge_rna, design)
fit_rna = glmFit(dge_rna, design)
rna_res = as.data.frame(glmLRT(fit_rna, coef="Age")) %>%
  merge(rinfo, ., by.x="Geneid", by.y=0) %>%
  mutate(class = factor(class, levels=re_order))

write.csv(rna_res, paste0(paths$out, "/supp/rna_seq_results.csv"), quote=F, row.names = F)

# Class boxes
fig$S3$rna_by_class = rna_res %>%
  left_join(stats_class, by="class") %>%
  dplyr::filter(n > 100) %>%
  mutate(selfish = ifelse(class %in% c("LINE", "SINE", "LTR", "DNA", "Retroposon"), "Selfish", "Non-selfish")) %>%
  ggplot(., aes(x=class, y=logFC, fill=selfish))+
  geom_boxplot(outlier.shape=NA)+
  geom_hline(yintercept=0, linetype="dashed", color="dark green")+
  scale_x_discrete(limits=rev)+
  scale_fill_manual(values=c("#1177dd", "#dd1111"))+
  coord_flip(ylim=c(-0.012, 0.012))+
  labs(y="Expression drift rate [logCPM/year]", fill="Repeat type")+
  theme(axis.title.y=element_blank())+
  ggtitle("Expression drift of RE classes")

# Bars
levs = stats_l1 %>% 
  dplyr::filter(n > 40) %>%
  arrange(mean_len) %>%
  pull(fam)
fig$S3$rna_l1 = rna_res %>%
  dplyr::filter(superf == "LINE/L1") %>%
  merge(stats_l1, ., by=c("fam")) %>%
  dplyr::filter(n>40)%>%
  mutate(fam = factor(fam, levels = levs)) %>%
  # mutate(fam = fct_reorder(fam, mean_len)) %>%
  mutate(padj = p.adjust(PValue, method="BH")) %>%
  ggplot(., aes(fam, logFC))+
  geom_bar(stat="identity")+
  scale_x_discrete(drop = FALSE)+
  theme(axis.text.x = element_text(angle=90, hjust=1),
        axis.title.x =  element_blank())+
  labs(y="Expression drift rate [logCPM/year]")
# # Smooth fit
# cpm2 = cpm(dge1, log = T)
# tmp = as.data.frame(t(cpm2[c("L1HS:L1:LINE", "L1PA2:L1:LINE", "L1PA3:L1:LINE", "L1PA4:L1:LINE"), ]))
# colnames(tmp) = c("L1HS", "L1PA2", "L1PA3", "L1PA4")
# all.equal(rownames(tmp), meta_rna$SRR_rna)
# p2 = cbind(meta_rna, tmp) %>%
#   pivot_longer(cols=c("L1HS", "L1PA2", "L1PA3", "L1PA4"), names_to="family", values_to = "cpm") %>%
#   ggplot(., aes(Age, cpm))+
#   geom_point(size=0.6)+
#   geom_smooth()+
#   geom_smooth(method="lm", se=F, color="red", linetype="dashed")+
#   facet_wrap(~family, scales = "free_y")+
#   stat_cor()
# p1+p2
# ggsave(paste0(paths$out, "/rna_seq_L1.png"), width=1.3*w, height=0.4*h, units="mm", dpi = 300)
# fig$S3$rna_l1 = p1

# Correlation between the two
rna_res2 = probe_info %>%
  dplyr::filter(!is.na(class)) %>%
  dplyr::select(class, superf, fam, age_coef, residuals) %>%
  group_by(class, superf, fam) %>%
  summarize(median_age_coef = median(age_coef), median_residual = median(residuals)) %>%
  merge(rna_res, ., by=c("fam"))
# Separate L1 and L2 because it makes a big difference
rna_res2$category = as.character(rna_res2$class.y)
rna_res2[rna_res2$class.y == "LINE", "category"] = as.character(rna_res2[rna_res2$class.y == "LINE", "superf.y"])
keep = table(rna_res2$category)
keep = setdiff(rownames(keep)[keep > 10], "Unknown")
keep = c("LINE/L1", "LINE/CR1", "SINE", "LTR", "DNA", "Satellite")
corP = function(x,y) cor.test(x,y)$p.value
# Correlation stats
stat_tmp = rna_res2 %>%
  dplyr::filter(category %in% keep) %>%
  group_by(category) %>%
  summarize(maxY = max(logFC), minX = min(median_age_coef),
            r = cor(median_age_coef, logFC), p = corP(median_age_coef, logFC)) %>%
  mutate(r = sprintf("r = %.3f", r)) %>%
  mutate(p = ifelse(p< 0.001, "p < 0.001", sprintf("p = %.3f", p))) %>%
  mutate(lab = paste(r, p, sep=", "))
fig$S3$met_rna_cor = rna_res2 %>%
  dplyr::filter(category %in% keep) %>%
  mutate(category = factor(category, levels = keep))%>%
  ggplot(., aes(median_age_coef, logFC))+
  geom_point(size=0.5)+
  geom_smooth(method="lm", se=F)+
  facet_wrap(~category, scales = "free", nrow=2)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  geom_text(data = stat_tmp, aes(x=minX*0.95, y = maxY*1.15, label=lab), hjust=0)+
  labs(x="Methylation drift rate [%/year]", y="Accessibility drift rate [logCPM/year]")

##### ATAC #####

load(paste0(paths$data, "/processed/atac_data.Rdata"))

# Fit EdgeR
all.equal(rownames(atac_re_tpm), rownames(rinfo_atac))
dge_atac = DGEList(atac_re_counts, samples = meta_atac)
design = model.matrix(~Age+Sex, data=meta_atac)
dge_atac = calcNormFactors(dge_atac)
dge_atac = estimateDisp(dge_atac, design)
fit_atac = glmFit(dge_atac, design)
atac_res = as.data.frame(glmLRT(fit_atac, coef="Age")) %>%
  merge(rinfo_atac, .,by=0)%>%
  mutate(class = factor(class, levels=re_order))

write.csv(atac_res, paste0(paths$out, "/supp/atac_seq_results.csv"), quote=F, row.names = F)

# Class boxes
fig$S3$atac_by_class = atac_res %>%
  left_join(stats_class, by="class") %>%
  dplyr::filter(n > 100) %>%
  mutate(selfish = ifelse(class %in% c("LINE", "SINE", "LTR", "DNA", "Retroposon"), "Selfish", "Non-selfish")) %>%
  ggplot(., aes(x=class, y=logFC, fill=selfish))+
  geom_boxplot(outlier.shape=NA)+
  geom_hline(yintercept=0, linetype="dashed", color="dark green")+
  scale_x_discrete(limits=rev)+
  scale_fill_manual(values=c("#1177dd", "#dd1111"))+
  coord_flip(ylim=c(-0.02, 0.02))+
  labs(y="Accessibility drift rate [logCPM/year]", fill="Repeat type")+
  theme(axis.title.y=element_blank())+
  ggtitle("Accessibility drift of RE classes")

# Bars
levs = stats_l1 %>% 
  dplyr::filter(n > 40) %>%
  arrange(mean_len) %>%
  pull(fam)
fig$S3$atac_l1 = atac_res %>%
  dplyr::filter(superf == "LINE/L1") %>%
  merge(stats_l1, ., by=c("superf", "fam")) %>%
  dplyr::filter(n>40)%>%
  mutate(fam = factor(fam, levels = levs)) %>%
  # mutate(fam = fct_reorder(fam, mean_len)) %>%
  mutate(padj = p.adjust(PValue, method="BH")) %>%
  ggplot(., aes(fam, logFC))+
  geom_bar(stat="identity")+
  scale_x_discrete(drop = FALSE)+
  theme(axis.text.x = element_text(angle=90, hjust=1),
        axis.title.x =  element_blank())+
  labs(y="Accessibility drift rate [logCPM/year]")

# Correlation between the two
atac_res2 = probe_info %>%
  dplyr::filter(!is.na(class)) %>%
  dplyr::select(class, superf, fam, age_coef, residuals) %>%
  group_by(class, superf, fam) %>%
  summarize(median_age_coef = median(age_coef), median_residual = median(residuals)) %>%
  merge(atac_res, ., by=c("class", "superf", "fam"))
# Separate L1 and L2 because it makes a big difference
atac_res2$category = as.character(atac_res2$class)
atac_res2[atac_res2$class == "LINE", "category"] = atac_res2[atac_res2$class == "LINE", "superf"]
keep = table(atac_res2$category)
keep = setdiff(rownames(keep)[keep > 10], "Unknown")
keep = c("LINE/L1", "LINE/L2", "LINE/CR1", 
         "SINE", "LTR", "DNA", "Satellite", "Simple_repeat", "tRNA", "snRNA")
corP = function(x,y) cor.test(x,y)$p.value
# Correlation stats
stat_tmp = atac_res2 %>%
  dplyr::filter(category %in% keep) %>%
  group_by(category) %>%
  summarize(maxY = max(logFC), minX = min(median_age_coef),
            r = cor(median_age_coef, logFC), p = corP(median_age_coef, logFC)) %>%
  mutate(r = sprintf("r = %.3f", r)) %>%
  mutate(p = ifelse(p< 0.001, "p < 0.001", sprintf("p = %.3f", p))) %>%
  mutate(lab = paste(r, p, sep=", "))

fig$S3$met_atac_cor = atac_res2 %>%
  dplyr::filter(category %in% keep) %>%
  mutate(category = factor(category, levels = keep))%>%
  ggplot(., aes(median_age_coef, logFC))+
  geom_point(size=0.5)+
  geom_smooth(method="lm", se=F)+
  facet_wrap(~category, scales = "free", nrow=2)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  geom_text(data = stat_tmp, aes(x=minX*0.95, y = maxY*1.15, label=lab), hjust=0)+
  labs(x="Methylation drift rate [%/year]", y="Accessibility drift rate [logCPM/year]")

##### Combine #####

atac_res_l1 = atac_res %>%
  dplyr::filter(superf == "LINE/L1") %>%
  merge(stats_l1, ., by=c("superf", "fam")) %>%
  dplyr::filter(n>40) %>%
  mutate(fam = factor(fam, levels = levs)) %>%
  # mutate(fam = fct_reorder(fam, mean_len)) %>%
  mutate(padj = p.adjust(PValue, method="BH"))
rna_res_l1 = rna_res %>%
  dplyr::filter(superf == "LINE/L1") %>%
  merge(stats_l1, ., by=c("superf", "fam")) %>%
  dplyr::filter(n>40) %>%
  mutate(fam = factor(fam, levels = levs)) %>%
  # mutate(fam = fct_reorder(fam, mean_len)) %>%
  mutate(padj = p.adjust(PValue, method="BH"))

fig$S3$l1_length = stats_l1 %>%
  dplyr::filter(n > 40) %>%
  mutate(fam = fct_reorder(fam, mean_len)) %>%
  ggplot(., aes(x=fam, y=mean_len, group=1))+
  geom_line()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+
  labs(y="Mean seq\nlength")
fig$S3$l1_length / fig$S3$rna_l1 / fig$S3$atac_l1 + plot_layout(heights = c(1,4,4))

#### PROBE LOCATIONS ####

##### Prepare consensus sequences #####

l1pa_cons = readDNAStringSet(paste0(paths$data, "/annotation/te_cons_repbase/l1pa_cons.fa"))

l1pa2_3p = readDNAStringSet(paste0(paths$data, "/annotation/te_cons_repbase/l1pa2_3p.fa"))
l1pa3_3p = readDNAStringSet(paste0(paths$data, "/annotation/te_cons_repbase/l1pa3_3p.fa"))
l1pa4_3p = readDNAStringSet(paste0(paths$data, "/annotation/te_cons_repbase/l1pa4_3p.fa"))

l1ma3_3p = readDNAStringSet(paste0(paths$data, "/annotation/te_cons_repbase/l1ma3_3p.fa"))
l1pa15_3p = readDNAStringSet(paste0(paths$data, "/annotation/te_cons_repbase/l1pa15_3p.fa"))
l1pa16_3p = readDNAStringSet(paste0(paths$data, "/annotation/te_cons_repbase/l1pa16_3p.fa"))
l1pb1_3p = readDNAStringSet(paste0(paths$data, "/annotation/te_cons_repbase/l1pb1_3p.fa"))

the1_int = readDNAStringSet(paste0(paths$data, "/annotation/te_cons_repbase/the1_int.fa"))
the1a_ltr = readDNAStringSet(paste0(paths$data, "/annotation/te_cons_repbase/the1a_ltr.fa"))
the1c_ltr = readDNAStringSet(paste0(paths$data, "/annotation/te_cons_repbase/the1c_ltr.fa"))

# Body and 3' should overlap for 150 bases from 5254 onward
substr(l1pa_cons, 5254, 5254+150)
substr(l1pa2_3p, 1, 150)
substr(l1pa3_3p, 1, 150)
substr(l1pa4_3p, 1, 150)
substr(l1ma3_3p, 1, 150)
substr(l1pa15_3p, 1, 150)
substr(l1pa16_3p, 1, 150)
substr(l1pb1_3p, 1, 150)
# Indeed it does, so lets join them
l1pa2_cons = DNAStringSet(paste0(substr(l1pa_cons, 1, 5253), l1pa2_3p))
l1pa3_cons = DNAStringSet(paste0(substr(l1pa_cons, 1, 5253), l1pa3_3p))
l1pa4_cons = DNAStringSet(paste0(substr(l1pa_cons, 1, 5253), l1pa4_3p))
l1ma3_cons = DNAStringSet(paste0(substr(l1pa_cons, 1, 5253), l1ma3_3p))
l1pa15_cons = DNAStringSet(paste0(substr(l1pa_cons, 1, 5253), l1pa15_3p))
l1pa16_cons = DNAStringSet(paste0(substr(l1pa_cons, 1, 5253), l1pa16_3p))
l1pb1_cons = DNAStringSet(paste0(substr(l1pa_cons, 1, 5253), l1pb1_3p))

# LTRs don't overlap for THE1A/C so just join
the1a_cons = DNAStringSet(paste0(the1a_ltr, the1_int))
the1c_cons = DNAStringSet(paste0(the1c_ltr, the1_int))

#Then load the ones that are good to go
l1hs_cons = readDNAStringSet(paste0(paths$data, "/annotation/te_cons_repbase/l1hs_cons.fa"))

##### Get probe pos on consensus #####

# Get some columns that were only in the original manifest by illumina
manifest = fread(paste0(paths$data, "/annotation/short_probe_info_450.csv"), data.table=F)
manifest = manifest %>%
  dplyr::filter(IlmnID %in% probe_info$IlmnID)
rownames(manifest) = manifest$IlmnID
manifest = manifest[probe_info$IlmnID, ]
probe_info$strand = ifelse(manifest$Strand == "F", "+", "-")
probe_info$type = manifest[probe_info$IlmnID, "Infinium_Design_Type"]
probe_info$alleleA_probe = manifest[probe_info$IlmnID, "AlleleA_ProbeSeq"]

find_on_cons = function(pinfo, cons) {
  # Define regions of interest and add some extra columns
  tmp = data.frame(
    chr = paste0("chr", pinfo$chr),
    start = ifelse(pinfo$strand == "+", pinfo$MAPINFO, pinfo$MAPINFO - 48),
    end = ifelse(pinfo$strand == "+", pinfo$MAPINFO + 49, pinfo$MAPINFO + 1),
    strand = pinfo$strand,
    type = pinfo$type,
    IlmnID = pinfo$IlmnID,
    cpgpos = pinfo$MAPINFO,
    fam = pinfo$fam,
    re_start = pinfo$re_start,
    re_end = pinfo$re_end,
    re_strand = pinfo$re_strand)
  tmp[tmp$type == "II" & tmp$strand == "+", "start"] = tmp[tmp$type == "II" & tmp$strand == "+", "start"] + 1
  tmp[tmp$type == "II" & tmp$strand == "+", "end"] = tmp[tmp$type == "II" & tmp$strand == "+", "end"] + 1
  tmp[tmp$type == "II" & tmp$strand == "-", "start"] = tmp[tmp$type == "II" & tmp$strand == "-", "start"] - 1
  tmp[tmp$type == "II" & tmp$strand == "-", "end"] = tmp[tmp$type == "II" & tmp$strand == "-", "end"] - 1
  # Get probe target sequences
  ranges = GRanges(
    seqnames = tmp$chr,
    ranges = IRanges(
      start = tmp$start, end = tmp$end),
    strand = tmp$strand)
  names(ranges) = tmp$IlmnID
  seqs = getSeq(BSgenome.Hsapiens.UCSC.hg19, ranges)
  als_probe_to_cons_fw = pairwiseAlignment(pattern = seqs[tmp$strand == tmp$re_strand], subject = cons, type="local")
  als_probe_to_cons_rv = pairwiseAlignment(pattern = reverseComplement(seqs[tmp$strand != tmp$re_strand]), subject = cons, type="local")
  tmp[tmp$strand == tmp$re_strand, "score_probe_to_con"] = als_probe_to_cons_fw@score
  tmp[tmp$strand != tmp$re_strand, "score_probe_to_con"] = als_probe_to_cons_rv@score
  # Get repeat sequences
  ranges = GRanges(
    seqnames = tmp$chr,
    ranges = IRanges(
      start = tmp$re_start, end = tmp$re_end),
    strand = tmp$re_strand)
  re_seqs = getSeq(BSgenome.Hsapiens.UCSC.hg19, ranges)
  # Probe to consensus alignment
  # als_re_to_cons = pairwiseAlignment(pattern = re_seqs, subject = cons, type="overlap")
  # tmp$score_re_to_cons = als_re_to_cons@score
  tmp[tmp$strand == tmp$re_strand, "pos_on_cons"] = start(subject(als_probe_to_cons_fw)) - start(Biostrings::pattern(als_probe_to_cons_fw))
  tmp[tmp$strand != tmp$re_strand, "pos_on_cons"] = start(subject(als_probe_to_cons_rv)) - start(Biostrings::pattern(als_probe_to_cons_rv)) + 50 # If rev cpg is at end of read, not start
  tmp$age_coef = pinfo$age_coef
  tmp$age_coef_adj = pinfo$residuals
  out = list(
    regions = tmp,
    probe_seqs = seqs,
    re_seqs = re_seqs,
    probe_to_cons_fw = als_probe_to_cons_fw,
    probe_to_cons_rv = als_probe_to_cons_rv) #re_to_cons = als_re_to_cons
  return(out)
}

# # Did some double checking over here
# probes_l1pa2 = probe_info %>%
#   dplyr::filter(fam == "L1PA2")
# test = find_on_cons(probes_l1pa2, DNAStringSet(l1pa2_cons))
# 
# # Objective: check if provided probe resembles extracted target seq
# ind=10
# test$probe_seqs[ind]
# c(test$probe_seqs[ind], reverseComplement(DNAStringSet(probes_l1pa2[ind, "alleleA_probe"])))
# 
# # Objective: find probe region in re_to_cons alignment
# ind = 2 # 2 has bad probe alignment
# sprintf("Probe strand: %s, RE strand: %s", test$regions[ind, "strand"], test$region[ind, "re_strand"])
# offset = test$regions[ind, "cpgpos"] - test$regions[ind, "re_start"]
# s1 = substr(test$re_seqs[ind], offset-48, offset+1)
# s2 = substr(test$probe_seqs[ind], 1, 50)
# c(DNAStringSet(s1), reverseComplement(DNAStringSet(s2)))
# # So now find this in the re_to_cons alignment
# test$re_to_cons[2]
# s3 = substr(pattern(test$re_to_cons)[2], offset-48, offset+1)
# s3
# # Copy and consensus are almost perfect matches so by using almost the same regions i should find the probe
# shift = 2
# s4 = substr(l1pa2_cons, offset-48+shift, offset+1+shift)
# c(DNAStringSet(s1), reverseComplement(DNAStringSet(s2)), s4)

probes_l1hs = subset(probe_info, fam == "L1HS")
probes_l1pa2 = subset(probe_info, fam == "L1PA2")
probes_l1pa3 = subset(probe_info, fam == "L1PA3")
probes_l1pa4 = subset(probe_info, fam == "L1PA4")
probes_the1a = subset(probe_info, fam == "THE1A")
probes_the1c = subset(probe_info, fam == "THE1C")
probes_l1ma3 = subset(probe_info, fam == "L1MA3")
probes_l1pa15 = subset(probe_info, fam == "L1PA15")
probes_l1pa16 = subset(probe_info, fam == "L1PA16")
probes_l1pb1 = subset(probe_info, fam == "L1PB1")

map_to_cons = list(
  l1hs = find_on_cons(probes_l1hs, l1hs_cons),
  l1pa2 = find_on_cons(probes_l1pa2, l1pa2_cons),
  l1pa3 = find_on_cons(probes_l1pa3, l1pa3_cons),
  l1pa4 = find_on_cons(probes_l1pa4, l1pa4_cons),
  l1ma3 = find_on_cons(probes_l1ma3, l1ma3_cons),
  l1pa15 = find_on_cons(probes_l1pa15, l1pa15_cons),
  l1pa16 = find_on_cons(probes_l1pa16, l1pa16_cons),
  l1pb1 = find_on_cons(probes_l1pb1, l1pb1_cons),
  the1a = find_on_cons(probes_the1a, the1a_cons),
  the1c = find_on_cons(probes_the1c, the1c_cons))

# Check mapping qualities, OK
ggplot(map_to_cons$l1hs$regions, aes(score_probe_to_con))+
  geom_histogram()
ggplot(map_to_cons$l1pa2$regions, aes(score_probe_to_con))+
  geom_histogram()
ggplot(map_to_cons$l1pa3$regions, aes(score_probe_to_con))+
  geom_histogram()
ggplot(map_to_cons$l1pa4$regions, aes(score_probe_to_con))+
  geom_histogram()

ggplot(map_to_cons$l1ma3$regions, aes(score_probe_to_con))+
  geom_histogram()
ggplot(map_to_cons$l1pa15$regions, aes(score_probe_to_con))+
  geom_histogram()
ggplot(map_to_cons$l1pa16$regions, aes(score_probe_to_con))+
  geom_histogram()
ggplot(map_to_cons$l1pb1$regions, aes(score_probe_to_con))+
  geom_histogram()

ggplot(map_to_cons$the1a$regions, aes(score_probe_to_con))+
  geom_histogram()
ggplot(map_to_cons$the1c$regions, aes(score_probe_to_con))+
  geom_histogram()

##### Make annotations plots #####

# ORF and UTR annotation from repbase embl files
l1hs_sections = data.frame(
  start = c(1, 908, 1988, 5813),
  end = c(907, 1921, 5812, width(l1hs_cons)),
  feature = c("5'UTR", "ORF1", "ORF2", "3'UTR"))
l1pa2_sections = data.frame(
  start = c(1, 1030, 2110, 5404),
  end = c(1029, 2046, 5403, width(l1pa2_cons)),
  feature = c("5'UTR", "ORF1", "ORF2", "3'UTR"))
l1pa3_sections = data.frame(
  start = c(1, 1030, 2110, 5404),
  end = c(1029, 2046, 5403, width(l1pa3_cons)),
  feature = c("5'UTR", "ORF1", "ORF2", "3'UTR"))
l1pa4_sections = data.frame(
  start = c(1, 1030, 2110, 5404),
  end = c(1029, 2046, 5403, width(l1pa4_cons)),
  feature = c("5'UTR", "ORF1", "ORF2", "3'UTR"))

l1ma3_sections = data.frame(
  start = c(1, 1030, 2110, 5404),
  end = c(1029, 2046, 5403, width(l1ma3_cons)),
  feature = c("5'UTR", "ORF1", "ORF2", "3'UTR"))
l1pa15_sections = data.frame(
  start = c(1, 1030, 2110, 5404),
  end = c(1029, 2046, 5403, width(l1pa15_cons)),
  feature = c("5'UTR", "ORF1", "ORF2", "3'UTR"))
l1pa16_sections = data.frame(
  start = c(1, 1030, 2110, 5404),
  end = c(1029, 2046, 5403, width(l1pa16_cons)),
  feature = c("5'UTR", "ORF1", "ORF2", "3'UTR"))
l1pb1_sections = data.frame(
  start = c(1, 1030, 2110, 5404),
  end = c(1029, 2046, 5403, width(l1pb1_cons)),
  feature = c("5'UTR", "ORF1", "ORF2", "3'UTR"))

the1a_sections = data.frame(
  start = c(1, width(the1a_ltr)+1),
  end = c(width(the1a_ltr), width(the1a_cons)),
  feature = c("LTR", "Int"))
the1c_sections = data.frame(
  start = c(1, width(the1c_ltr)+1),
  end = c(width(the1c_ltr), width(the1c_cons)),
  feature = c("LTR", "Int"))
# Middle of sections for text
l1hs_sections$mid = (l1hs_sections$end + l1hs_sections$start) / 2
l1pa2_sections$mid = (l1pa2_sections$end + l1pa2_sections$start) / 2
l1pa3_sections$mid = (l1pa3_sections$end + l1pa3_sections$start) / 2
l1pa4_sections$mid = (l1pa4_sections$end + l1pa4_sections$start) / 2

l1ma3_sections$mid = (l1ma3_sections$end + l1ma3_sections$start) / 2
l1pa15_sections$mid = (l1pa15_sections$end + l1pa15_sections$start) / 2
l1pa16_sections$mid = (l1pa16_sections$end + l1pa16_sections$start) / 2
l1pb1_sections$mid = (l1pb1_sections$end + l1pb1_sections$start) / 2

the1a_sections$mid = (the1a_sections$end + the1a_sections$start) / 2
the1c_sections$mid = (the1c_sections$end + the1c_sections$start) / 2
# Plots
l1hs_bar = ggplot(l1hs_sections, aes(x = mid, y=0, xmin=start, xmax=end, ymin=-1, ymax=1, fill=feature, label=feature))+
  geom_rect()+
  geom_text()+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())+
  labs(x="Consensus position")
l1pa2_bar = ggplot(l1pa2_sections, aes(x = mid, y=0, xmin=start, xmax=end, ymin=-1, ymax=1, fill=feature, label=feature))+
  geom_rect()+
  geom_text()+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())+
  labs(x="Consensus position")
l1pa3_bar = ggplot(l1pa3_sections, aes(x = mid, y=0, xmin=start, xmax=end, ymin=-1, ymax=1, fill=feature, label=feature))+
  geom_rect()+
  geom_text()+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())+
  labs(x="Consensus position")
l1pa4_bar = ggplot(l1pa4_sections, aes(x = mid, y=0, xmin=start, xmax=end, ymin=-1, ymax=1, fill=feature, label=feature))+
  geom_rect()+
  geom_text()+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())+
  labs(x="Consensus position")

l1ma3_bar = ggplot(l1ma3_sections, aes(x = mid, y=0, xmin=start, xmax=end, ymin=-1, ymax=1, fill=feature, label=feature))+
  geom_rect()+
  geom_text()+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())+
  labs(x="Consensus position")
l1pa15_bar = ggplot(l1pa15_sections, aes(x = mid, y=0, xmin=start, xmax=end, ymin=-1, ymax=1, fill=feature, label=feature))+
  geom_rect()+
  geom_text()+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())+
  labs(x="Consensus position")
l1pa16_bar = ggplot(l1pa16_sections, aes(x = mid, y=0, xmin=start, xmax=end, ymin=-1, ymax=1, fill=feature, label=feature))+
  geom_rect()+
  geom_text()+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())+
  labs(x="Consensus position")
l1pb1_bar = ggplot(l1pb1_sections, aes(x = mid, y=0, xmin=start, xmax=end, ymin=-1, ymax=1, fill=feature, label=feature))+
  geom_rect()+
  geom_text()+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())+
  labs(x="Consensus position")


the1a_bar = ggplot(the1a_sections, aes(x = mid, y=0, xmin=start, xmax=end, ymin=-1, ymax=1, fill=feature, label=feature))+
  geom_rect()+
  geom_text()+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())+
  labs(x="Consensus position")
the1c_bar = ggplot(the1c_sections, aes(x = mid, y=0, xmin=start, xmax=end, ymin=-1, ymax=1, fill=feature, label=feature))+
  geom_rect()+
  geom_text()+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())+
  labs(x="Consensus position")

##### Final plots #####

fig$S1$cons_loc_l1hs_top = ggplot(map_to_cons$l1hs$regions, aes(pos_on_cons, age_coef))+
  geom_point(size=1)+
  xlim(0, width(l1hs_cons))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y = "Methylation drift rate [%/year]")+
  ggtitle("L1HS")
fig$S1$cons_loc_l1hs_bot = l1hs_bar
fig$S1$cons_loc_l1hs = fig$S1$cons_loc_l1hs_top / l1hs_bar + plot_layout(heights = c(4,1))

fig$S1$cons_loc_l1pa2_top = ggplot(map_to_cons$l1pa2$regions, aes(pos_on_cons, age_coef))+
  geom_point(size=1)+
  xlim(0, width(l1pa2_cons))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y = "Methylation drift rate [%/year]")+
  ggtitle("L1PA2")
fig$S1$cons_loc_l1pa2_bot = l1pa2_bar
fig$S1$cons_loc_l1pa2 = fig$S1$cons_loc_l1pa2_top / l1pa2_bar + plot_layout(heights = c(4,1))

fig$S1$cons_loc_l1pa3_top = ggplot(map_to_cons$l1pa3$regions, aes(pos_on_cons, age_coef))+
  geom_point(size=1)+
  xlim(0, width(l1pa3_cons))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y = "Methylation drift rate [%/year]")+
  ggtitle("L1PA3")
fig$S1$cons_loc_l1pa3_bot = l1pa3_bar
fig$S1$cons_loc_l1pa3 = fig$S1$cons_loc_l1pa3_top / l1pa3_bar + plot_layout(heights = c(4,1))

fig$S1$cons_loc_l1pa4_top = ggplot(map_to_cons$l1pa4$regions, aes(pos_on_cons, age_coef))+
  geom_point(size=1)+
  xlim(0, width(l1pa4_cons))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y = "Methylation drift rate [%/year]")+
  ggtitle("L1PA4")
fig$S1$cons_loc_l1pa4_bot = l1pa4_bar
fig$S1$cons_loc_l1pa4 = fig$S1$cons_loc_l1pa4_top / l1pa4_bar + plot_layout(heights = c(4,1))

fig$S1$cons_loc_l1ma3_top = ggplot(map_to_cons$l1ma3$regions, aes(pos_on_cons, age_coef))+
  geom_point(size=1)+
  xlim(0, width(l1ma3_cons))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y = "Methylation drift rate [%/year]")+
  ggtitle("L1MA3")
fig$S1$cons_loc_l1ma3_bot = l1ma3_bar
fig$S1$cons_loc_l1ma3 = fig$S1$cons_loc_l1ma3_top / l1ma3_bar + plot_layout(heights = c(4,1))

fig$S1$cons_loc_l1pa15_top = ggplot(map_to_cons$l1pa15$regions, aes(pos_on_cons, age_coef))+
  geom_point(size=1)+
  xlim(0, width(l1pa15_cons))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y = "Methylation drift rate [%/year]")+
  ggtitle("L1PA15")
fig$S1$cons_loc_l1pa15_bot = l1pa15_bar
fig$S1$cons_loc_l1pa15 = fig$S1$cons_loc_l1pa15_top / l1pa15_bar + plot_layout(heights = c(4,1))

fig$S1$cons_loc_l1pa16_top = ggplot(map_to_cons$l1pa16$regions, aes(pos_on_cons, age_coef))+
  geom_point(size=1)+
  xlim(0, width(l1pa16_cons))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y = "Methylation drift rate [%/year]")+
  ggtitle("L1PA16")
fig$S1$cons_loc_l1pa16_bot = l1pa16_bar
fig$S1$cons_loc_l1pa16 = fig$S1$cons_loc_l1pa16_top / l1pa16_bar + plot_layout(heights = c(4,1))

fig$S1$cons_loc_l1pb1_top = ggplot(map_to_cons$l1pb1$regions, aes(pos_on_cons, age_coef))+
  geom_point(size=1)+
  xlim(0, width(l1pb1_cons))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y = "Methylation drift rate [%/year]")+
  ggtitle("L1PB1")
fig$S1$cons_loc_l1pb1_bot = l1pb1_bar
fig$S1$cons_loc_l1pb1 = fig$S1$cons_loc_l1pb1_top / l1pb1_bar + plot_layout(heights = c(4,1))

fig$S1$cons_loc_the1a_top = ggplot(map_to_cons$the1a$regions, aes(pos_on_cons, age_coef))+
  geom_point(size=1)+
  xlim(0, width(the1a_cons))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y = "Methylation drift rate [%/year]")+
  ggtitle("THE1A")
fig$S1$cons_loc_the1a_bot = the1a_bar
fig$S1$cons_loc_the1a = fig$S1$cons_loc_the1a_top / the1a_bar + plot_layout(heights = c(4,1))

fig$S1$cons_loc_the1c_top = ggplot(map_to_cons$the1c$regions, aes(pos_on_cons, age_coef))+
  geom_point(size=1)+
  xlim(0, width(the1c_cons))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y = "Methylation drift rate [%/year]")+
  ggtitle("THE1C")
fig$S1$cons_loc_the1c_bot = the1c_bar
fig$S1$cons_loc_the1c = fig$S1$cons_loc_the1c_top / the1c_bar + plot_layout(heights = c(4,1))

fig$S1$cons_loc_assembled = wrap_plots(fig$S1[seq(5,length(fig$S1),3)], ncol=2)

#### CHECKPOINT ####

# save.image(paste0(paths$out, "/checkpoint_p1.Rdata"))
load(paste0(paths$out, "/checkpoint_p1.Rdata"))

#### MOTIF ENRICHMENT ####

##### Exploration #####

l1_young = c("L1HS", "L1PA2", "L1PA3", "L1PA4")
l1_older = c("L1PB1", "L1PA15", "L1PA16", "L1M1", "L1M2", "L1MA3")
l1_oi = c(l1_young, l1_older)
probes_l1oi = probe_info %>%
  dplyr::filter(fam %in% l1_oi)
table(probes_l1oi$fam)

ggplot(probes_l1oi, aes(fam, age_coef))+
  geom_boxplot(aes(fill=after_stat(middle)))+
  scale_fill_gradient2(low="blue", high="red")
ggplot(probes_l1oi, aes(fam, residuals))+
  geom_boxplot(aes(fill=after_stat(middle)))+
  scale_fill_gradient2(low="blue", high="red")

ggplot(probes_l1oi, aes(age_coef))+
  geom_histogram(bins=100)
ggplot(probes_l1oi, aes(residuals))+
  geom_histogram(bins=100)

##### Get fixed size flank #####

flankn = 250
ranges = GRanges(
  seqnames = paste0("chr", probes_l1oi$chr),
  ranges = IRanges(start = probes_l1oi$MAPINFO - flankn, end = probes_l1oi$MAPINFO + flankn),
  strand = probes_l1oi$strand)
names(ranges) = probes_l1oi$IlmnID
flank_seqs = getSeq(BSgenome.Hsapiens.UCSC.hg19, ranges)
all.equal(names(flank_seqs), probes_l1oi$IlmnID)
probes_l1oi$flank_seq = as.character(flank_seqs)

##### Write young vs old seqs #####

write.fasta(sequences = as.list(probes_l1oi[probes_l1oi$fam %in% l1_older, "flank_seq"]), 
            names = probes_l1oi[probes_l1oi$fam %in% l1_older, "IlmnID"],
            file = paste0(paths$out, "/motif/flank/flank_seq_old.fa"))
write.fasta(sequences = as.list(probes_l1oi[probes_l1oi$fam %in% l1_young, "flank_seq"]), 
            names = probes_l1oi[probes_l1oi$fam %in% l1_young, "IlmnID"],
            file = paste0(paths$out, "/motif/flank/flank_seq_young.fa"))

##### Write within family comparison seqs #####

for (fm in l1_young) {
  top50 = probes_l1oi %>% # top drift in negative sense, so min not max
    dplyr::filter(fam == fm) %>%
    slice_min(order_by=residuals, prop=0.5)
  bot50 = probes_l1oi %>%
    dplyr::filter(fam == fm) %>%
    slice_max(order_by=residuals, prop=-0.5)
  stopifnot(length(intersect(top50$IlmnID, bot50$IlmnID)) == 0)
  write.fasta(sequences = as.list(top50$flank_seq), 
              names = top50$IlmnID,
              file = sprintf("%s/motif/flank/flank_seq_%s_top50.fa", paths$out, fm))
  write.fasta(sequences = as.list(bot50$flank_seq), 
              names = bot50$IlmnID,
              file = sprintf("%s/motif/flank/flank_seq_%s_bot50.fa", paths$out, fm))
}

##### Read SEA results #####

sea = list()
sea_sites = list()
sea_seqs = list()

# Only L1PA2 and L1PA4 have any enriched motif in top50 vs bot50
sea$y_vs_o = fread(cmd=sprintf("grep -v '^#' %s", paste0(paths$out, "/motif/results/sea_y_vs_o_jaspar.tsv")), data.table = F)
sea$l1pa2 = fread(cmd=sprintf("grep -v '^#' %s", paste0(paths$out, "/motif/results/sea_l1pa2.tsv")), data.table = F)
sea$l1pa4 = fread(cmd=sprintf("grep -v '^#' %s", paste0(paths$out, "/motif/results/sea_l1pa4.tsv")), data.table = F)

sea_sites$y_vs_o = fread(cmd=sprintf("grep -v '^#' %s", paste0(paths$out, "/motif/sites/y_vs_o_jaspar.tsv")), data.table = F)
sea_sites$l1pa2 = fread(cmd=sprintf("grep -v '^#' %s", paste0(paths$out, "/motif/sites/l1pa2.tsv")), data.table = F)
sea_sites$l1pa4 = fread(cmd=sprintf("grep -v '^#' %s", paste0(paths$out, "/motif/sites/l1pa4.tsv")), data.table = F)

sea_seqs$y_vs_o = fread(cmd=sprintf("grep -v '^#' %s", paste0(paths$out, "/motif/seqs/y_vs_o_jaspar.tsv")), data.table = F)
sea_seqs$l1pa2 = fread(cmd=sprintf("grep -v '^#' %s", paste0(paths$out, "/motif/seqs/l1pa2.tsv")), data.table = F)
sea_seqs$l1pa4 = fread(cmd=sprintf("grep -v '^#' %s", paste0(paths$out, "/motif/seqs/l1pa4.tsv")), data.table = F)

for (comp in names(sea)) {
  sea[[comp]]$uniqID = sea[[comp]]$ALT_ID
  sea[[comp]][sea[[comp]]$uniqID == "", "uniqID"] = sea[[comp]][sea[[comp]]$uniqID == "", "ID"]
  sea[[comp]]$uniqID = toupper(sea[[comp]]$uniqID)
  sea[[comp]]$uniqID = str_remove(sea[[comp]]$uniqID, "_DBD.*")
  sea[[comp]]$uniqID = str_remove(sea[[comp]]$uniqID, "_FULL.*")
  sea[[comp]]$uniqID = str_remove(sea[[comp]]$uniqID, "_PRIMARY.*")
  sea[[comp]]$uniqID = str_remove(sea[[comp]]$uniqID, "_SECONDARY.*")
  sea[[comp]]$uniqID = str_remove(sea[[comp]]$uniqID, "_FULL.*")
  sea[[comp]]$uniqID = str_remove(sea[[comp]]$uniqID, "_MOUSE.*")
}

length(unique(sea$y_vs_o$uniqID))

##### Filter very high in fg, very low in bg #####

# Restrict to motifs found in almost all fg, almost no bg
colnames(sea$y_vs_o) = str_replace(colnames(sea$y_vs_o), "%", "_perc")
ggplot(sea$y_vs_o, aes(x=TP_perc, y=FP_perc))+
  geom_density2d()
ggplot(sea$y_vs_o, aes(x=TP_perc))+
  geom_histogram()
ggplot(sea$y_vs_o, aes(x=FP_perc))+
  geom_histogram()
View(sea$y_vs_o)
sea_res2 = sea$y_vs_o %>%
  dplyr::filter(TP_perc > 50) %>% # Was 50
  dplyr::filter(FP_perc < 20) # Was 20
# test = sea_res %>%
#   dplyr::filter(!ID %in% sea_res2$ID)
sum(probes_l1oi$fam %in% l1_young)
sum(probes_l1oi$fam %in% l1_older)
length(unique(sea_res2$uniqID))

##### Filter based on effect within young L1s #####

sea_seqs2 = sea_seqs$y_vs_o %>%
  dplyr::filter(motif_ID %in% sea_res2$ID)
probes_l1oi$young = probes_l1oi$fam %in% l1_young

motif_ids = unique(sea_seqs2$motif_ID)
res = data.frame()
# i=1047
for (i in 1:length(motif_ids)) {
  positive_seqs = sea_seqs2 %>%
    dplyr::filter(motif_ID == motif_ids[i]) %>%
    pull(seq_ID)
  tmp = probes_l1oi %>%
    mutate(matched = IlmnID %in% positive_seqs)
  tmp$group = "Old"
  tmp[tmp$young & tmp$matched, "group"] = "Young_matched"
  tmp[tmp$young & !tmp$matched, "group"] = "Young_not_matched"
  res[i, "n_old"] = sum(tmp$group == "Old")
  res[i, "n_young_matched"] = sum(tmp$group == "Young_matched")
  res[i, "n_young_unmatched"] = sum(tmp == "Young_not_matched")
  res[i, "motif_ID"] = motif_ids[i]
  res[i, "medianCoef_Old"] = median(tmp[!tmp$young, "age_coef"])
  res[i, "medianCoef_Young"] = median(tmp[tmp$young, "age_coef"])
  res[i, "medianCoef_Young_matched"] = median(tmp[tmp$group == "Young_matched", "age_coef"])
  res[i, "medianCoef_Young_not_matched"] = median(tmp[tmp$group == "Young_not_matched", "age_coef"])
  res[i, "medianAdjCoef_Old"] = median(tmp[!tmp$young, "residuals"])
  res[i, "medianAdjCoef_Young"] = median(tmp[tmp$young, "residuals"])
  res[i, "medianAdjCoef_Young_matched"] = median(tmp[tmp$group == "Young_matched", "residuals"])
  res[i, "medianAdjCoef_Young_not_matched"] = median(tmp[tmp$group == "Young_not_matched", "residuals"])
  # Young matched vs young unmatched
  tmp1 = dplyr::filter(tmp, young)
  if (length(unique(tmp1$group)) > 1) {
    # Using lm, accounting for family
    test = lm(age_coef ~ 0+fam+matched, data=tmp1)
    res[i, "effect_coef_matched_in_young"] = test$coefficients["matchedTRUE"]
    res[i, "pval_coef_matched_in_young"] = summary(test)$coefficients[,4]["matchedTRUE"]
    test = lm(residuals ~ 0+fam+matched, data=tmp1)
    res[i, "effect_adjCoef_matched_in_young"] = test$coefficients["matchedTRUE"]
    res[i, "pval_adjCoef_matched_in_young"] = summary(test)$coefficients[,4]["matchedTRUE"]
    # Using wilcox (not accounting for family)
    test = wilcox.test(age_coef ~ matched, data=tmp1)
    res[i, "pval_coef_matched_in_young_wilcox"] = test$p.value
    test = wilcox.test(residuals ~ matched, data=tmp1)
    res[i, "pval_adjCoef_matched_in_young_wilcox"] = test$p.value
  }
}
res = res %>%
  mutate(padj_coef_matched_in_young = p.adjust(pval_coef_matched_in_young, method="BH"), .after="pval_coef_matched_in_young")%>%
  mutate(padj_adjCoef_matched_in_young = p.adjust(pval_adjCoef_matched_in_young, method="BH"), .after="pval_adjCoef_matched_in_young") %>%
  mutate(padj_coef_matched_in_young_wilcox = p.adjust(pval_coef_matched_in_young_wilcox, method="BH"), .after="pval_coef_matched_in_young_wilcox")%>%
  mutate(padj_adjCoef_matched_in_young_wilcox = p.adjust(pval_adjCoef_matched_in_young_wilcox, method="BH"), .after="pval_adjCoef_matched_in_young_wilcox")
res = sea_res2 %>%
  merge(res, ., by.x="motif_ID", by.y="ID") %>%
  mutate(MvU = medianCoef_Young_matched - medianCoef_Young_not_matched) %>%
  mutate(MvU_adj = medianAdjCoef_Young_matched - medianAdjCoef_Young_not_matched)

# Compare pvalues obtained with different methods
ggplot(res, aes(-log10(pval_coef_matched_in_young), -log10(pval_adjCoef_matched_in_young)))+
  geom_point()
ggplot(res, aes(-log10(pval_coef_matched_in_young), -log10(pval_coef_matched_in_young_wilcox)))+
  geom_point()
ggplot(res, aes(-log10(pval_adjCoef_matched_in_young), -log10(pval_adjCoef_matched_in_young_wilcox)))+
  geom_point()
# Compare effect estimates obtained with different methods
ggplot(res, aes(MvU, MvU_adj))+
  geom_point()
ggplot(res, aes(effect_coef_matched_in_young, effect_adjCoef_matched_in_young))+
  geom_point()
ggplot(res, aes(MvU, effect_coef_matched_in_young))+
  geom_point()
ggplot(res, aes(MvU_adj, effect_adjCoef_matched_in_young))+
  geom_point()

# Investigate thresholds
res2 = res %>%
  dplyr::filter(pval_coef_matched_in_young_wilcox < 0.05) %>%
  dplyr::filter(pval_adjCoef_matched_in_young_wilcox < 0.05) %>%
  mutate(isFox = grepl("fox", motif_ID, ignore.case = T) | grepl("fox", ALT_ID, ignore.case = T)) %>%
  arrange(MvU_adj)
ggplot(res2, aes(effect_coef_matched_in_young, effect_adjCoef_matched_in_young))+
  geom_point()
ggplot(res2, aes(pval_coef_matched_in_young_wilcox, pval_adjCoef_matched_in_young_wilcox))+
  geom_point()
ggplot(res2, aes(MvU_adj, -log(pval_adjCoef_matched_in_young_wilcox)))+
  geom_point()
 
res_uniq = distinct(res2, uniqID, .keep_all = T)
# Split repressive and activating
res_repr = res_uniq %>%
  dplyr::filter(MvU_adj > 0)
res_act = res_uniq %>%
  dplyr::filter(MvU_adj < 0)
res_nofox = res_repr[!res_repr$isFox, ]

plot_motif_rate_association = function(probes_l1, sea_seqs, res2, tf_name) {
  tf_id = res2[res2$uniqID == tf_name, "motif_ID"]
  probes_l1 %>%
    mutate(matched = IlmnID %in% sea_seqs[sea_seqs$motif_ID == tf_id, "seq_ID"]) %>% #nrow(res2)
    mutate(Group = factor(paste0(young, "_", matched))) %>%
    mutate(Group = fct_recode(Group, "OldL1" = "FALSE_TRUE", "OldL1" = "FALSE_FALSE",
                              "YoungL1_with_motif" = "TRUE_TRUE", "YoungL1_without_motif" = "TRUE_FALSE")) %>%
    ggplot(., aes(Group, residuals, fill=Group))+
    geom_boxplot(outlier.shape = NA)+
    coord_cartesian(ylim=c(-0.07, 0.05))+
    scale_fill_manual(values=c("#6C5E6B", "#F65653", "#00729E"))+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5))+
    labs(y="Adj. methylation drift rate [%/yr]")+
    ggtitle(tf_name)
}
fig$S4$AR = plot_motif_rate_association(probes_l1oi, sea_seqs$y_vs_o, res_uniq, "AR")+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
fig$S4$ARNTL = plot_motif_rate_association(probes_l1oi, sea_seqs$y_vs_o, res_uniq, "ARNTL")+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
fig$S4$FOXO1 = plot_motif_rate_association(probes_l1oi, sea_seqs$y_vs_o, res_uniq, "FOXO1")+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
fig$S4$HIC2 = plot_motif_rate_association(probes_l1oi, sea_seqs$y_vs_o, res_uniq, "HIC2")+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
fig$S4$NFKB2 = plot_motif_rate_association(probes_l1oi, sea_seqs$y_vs_o, res_uniq, "NFKB2")+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
fig$S4$NR1I2 = plot_motif_rate_association(probes_l1oi, sea_seqs$y_vs_o, res_uniq, "NR1I2")+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
fig$S4$RXRA = plot_motif_rate_association(probes_l1oi, sea_seqs$y_vs_o, res_uniq, "RXRA::VDR")+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
fig$S4$SP5 = plot_motif_rate_association(probes_l1oi, sea_seqs$y_vs_o, res_uniq, "SP5")+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

fig$S4$AR + fig$S4$ARNTL + fig$S4$FOXO1 + fig$S4$HIC2 + fig$S4$NFKB2 + fig$S4$NR1I2 + fig$S4$RXRA + fig$S4$SP5 + plot_layout(guides="collect")
ggsave(paste0(paths$out, "/motif_association.png"), width=1*w, height=0.4*h, units="mm", dpi = 300)
write.csv(res_uniq, paste0(paths$out, "/motif_association.csv"))

##### Motif analysis within young L1s #####

# Figures assembled by hand from SEA results

#### SAVE ####

write.csv(probe_info, paste0(paths$out, "/supp/probe_info.csv"), row.names = F, quote=F)
save(fig, file=paste0(paths$out, "/figures.Rdata"))
