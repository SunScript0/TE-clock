library(caret)
library(glmnet)
library(rjson)
library(riverplot)
library(patchwork)
library(sp)
library(MASS)

#### DONE ####

##### Utility #####

logit = function(data, pseudoc) {
  return(log(data+pseudoc) - log(1-data+pseudoc))
}

split_rep_id = function(rep_id) {
  # Split rep_id into superf, fam, id properly
  tmp = rep_id
  rep_id = str_extract(tmp, "[[:digit:]]+$")
  tmp = str_replace(tmp, "/[[:digit:]]+$", "")
  fam = str_extract(tmp, "[^/]+$")
  tmp = str_replace(tmp, "/[^/]+$", "")
  superf = tmp
  out = data.frame(
    superf = superf,
    fam = fam,
    rep_id = rep_id
  )
  return(out)
}

plot_nfeat_vs_mincov = function(cov) {
  ths = c(seq(0, 19, 1), seq(20, 100, 10))
  res = data.frame(
    ID = character(),
    th = numeric(),
    passing = numeric()
  )
  for (th in ths) {
    passing = rowSums(cov >= th)
    ids = names(passing)
    tmp2 = data.frame(
      ID = ids,
      th = rep(th, length(ids)),
      passing = passing
    )
    res = rbind(res, tmp2)
  }
  med = res %>%
    group_by(th) %>%
    summarize(mean = mean(passing))
  res[res$passing == 0, "passing"] = NA
  p = ggplot(res, aes(x=th, y=passing))+
    geom_line(color="gray", aes(group=ID))+
    geom_line(data=med, aes(x=th, y=mean), color="red")+
    scale_y_log10()+
    labs(x="Minimum coverage", y="Features passing")
  return(p)
}

plot_cors_vs_mincov = function(data) {
  cov_mean = colMeans(data$cov)
  ths = seq(0, 100, 5)
  res = data.frame(
    th = numeric(),
    cov_cors = numeric(),
    met_cors = numeric()
  )
  for (th in ths) {
    cov_cor = cor(t(data$cov[, cov_mean > th]), use = "pairwise.complete.obs", method = "spearman")
    cov_cor = cov_cor[upper.tri(cov_cor, diag=F)]
    tmp2 = data$met[, cov_mean > th] / data$cov[, cov_mean > th]
    met_cor = cor(t(tmp2), use = "pairwise.complete.obs", method = "spearman")
    met_cor = met_cor[upper.tri(met_cor, diag=F)]
    res = rbind(res, data.frame(th = rep(th, length(cov_cor)), cov_cors = cov_cor, met_cors = met_cor))
  }
  p1 = ggplot(res, aes(x=th, y=cov_cors, group=th))+
    geom_boxplot()+
    labs(x="Minimum coverage", y="Correlation of coverage\n between samples")
  p2 = ggplot(res, aes(x=th, y=met_cors, group=th))+
    geom_boxplot()+
    labs(x="Minimum coverage", y="Correlation of methylation\n between samples")
  return(list(p1, p2))
}

outlier_detection = function(data, meta, group = "study", subgroups = NULL,
                             savepath = NULL, w_in = 3.1, h_in = 2.4) {
  if (any(rownames(data) != rownames(meta))) {
    stop("data columns do not match meta rows")
  }
  meta$outlier = F
  ps = list()
  for (g in unique(meta[, group])) { # Main group (e.g. study)
    sset_meta = meta[meta[group] == g, ]
    sset_data = data[rownames(sset_meta), ]
    pca = prcomp(sset_data, center=T, scale.=T)
    pca = merge(pca$x[, c("PC1", "PC2")], sset_meta, by=0)
    if (is.null(subgroups)) {
      pca$subgroup = "none"
    } else {
      pca$subgroup = do.call(paste, sset_meta[, subgroups]) 
    }
    for (sg in unique(pca$subgroup)) { # Secondary group
      sset_pca = pca[pca$subgroup == sg, ]
      p = ggplot(sset_pca, aes(x=PC1, y=PC2))+
        geom_point()+
        stat_ellipse(level=0.95)
      ellipse = ggplot_build(p)$data[[2]]
      if (any(is.na(ellipse$x))) {
        print(paste("No ellipse: all samples in", g, sg, "pass"))
        outs = c()
      } else {
        outs = sset_pca[!point.in.polygon(sset_pca$PC1, sset_pca$PC2, ellipse$x, ellipse$y), "Row.names"]
      }
      pca[pca$Row.names %in% outs, "outlier"] = T
      meta[outs, "outlier"] = T
    }
    ps[[g]] = ggplot(pca, aes(x=PC1, y=PC2, color=outlier))+
      geom_point()+
      ggtitle(g)
  }
  if (!is.null(savepath)) {
    ps = align_patches(ps)
    pdf(savepath, width=w_in, height=h_in)
    print(ps)
    dev.off()
  }
  return(meta$outlier)
}

##### Clocks #####

rmse = function(real, pred) {sqrt(mean((pred - real)^2))}
mae = function(real, pred) {median(abs(pred - real))}

ageT = function(age, adult_age) {
  out = (age - adult_age) / (adult_age + 1)
  out[age <= adult_age] = (log(age+1)-log(adult_age+1))[age <= adult_age]
  return(out)
}
ageTinv = function(t_age, adult_age) {
  out = adult_age*t_age + t_age + adult_age
  out[!is.na(t_age) & t_age <= 0] = (exp(t_age)*(adult_age+1)-1)[!is.na(t_age) & t_age <= 0]
  return(out)
}

make_groups = function(ngroups, labels) {
  nsamp = length(labels)
  labels_sorted_ind = order(labels)
  groups = c()
  for (i in 1:nsamp){
    s_ind = labels_sorted_ind[i]
    groups[s_ind] = i %% ngroups
  }
  return(groups+1)
}

ncv = function(data, labels, seed=sample(1:1000, 1), k_outer=10,
               cv_function, cv_args = list()) {
  set.seed(seed)
  # Fold indices
  fold_ind = createFolds(labels, k=k_outer, list=F)
  # Results table
  res = data.frame()
  for (i in 1:k_outer) {
    # Train-test partition
    train_data = data[fold_ind != i, ]
    train_labels = labels[fold_ind != i]
    test_data = data[fold_ind == i, ]
    test_labels = labels[fold_ind == i]
    # Inner CV
    cv_fit = cv_function(train_data, train_labels, cv_args)
    preds = predict(cv_fit, test_data)
    tmp = cbind(test_labels, preds, rep(i, length(preds)))
    rownames(tmp) = rownames(test_data)
    res = rbind(res, tmp)
  }
  final_fit = cv_function(data, labels, cv_args)
  colnames(res) = c("age", "preds", "fold")
  res$fold = as.factor(res$fold)
  out = list(
    preds = res,
    fit = final_fit,
    seed = seed
  )
  return(out)
}

ncv_group = function(data, labels, groups, seed=sample(1:1000, 1),
                     cv_function, cv_args = list()) {
  set.seed(seed)
  # Results table
  res = data.frame()
  for (i in 1:length(unique(groups))) {
    # Train-test partition
    train_data = data[groups != i, ]
    train_labels = labels[groups != i]
    train_groups = groups[groups != i]
    test_data = data[groups == i, ]
    test_labels = labels[groups == i]
    # Inner CV
    cv_fit = cv_function(train_data, train_labels, train_groups, cv_args)
    preds = predict(cv_fit, test_data)
    tmp = cbind(test_labels, preds, rep(i, length(preds)))
    rownames(tmp) = rownames(test_data)
    res = rbind(res, tmp)
  }
  final_fit = cv_function(data, labels, groups, cv_args)
  colnames(res) = c("age", "preds", "fold")
  res$fold = as.factor(res$fold)
  out = list(
    preds = res,
    fit = final_fit,
    seed = seed
  )
  return(out)
}

simple_split = function(data, labels, seed=sample(1:1000, 1),
                        cv_function, cv_args = list(), 
                        filt_function=NULL, filt_args=list(), 
                        clust_function=NULL, clust_args=list()) {
  set.seed(seed)
  # Train-test partition
  train_ind = createDataPartition(labels, p = 0.8, list = FALSE)
  train_data = data[train_ind, ]
  train_labels = labels[train_ind]
  test_data = data[-train_ind, ]
  test_labels = labels[-train_ind]
  cat(sprintf("Partitioned data: train(%i, %i), test(%i, %i)", 
              nrow(train_data), ncol(train_data), nrow(test_data), ncol(test_data)), "\n")
  # # Cluster
  # if (!is.null(clust_function)) {
  #   re_info$clust = clust_function(train_data, re_info, clust_args)
  #   train_data=collapse_to_clusters(train_data, re_info)$data
  #   test_data=collapse_to_clusters(test_data, re_info)$data
  #   re_info = re_info %>%
  #     group_by(clust) %>%
  #     summarize(n=sum(n))
  #   cat(sprintf("Clustered features into %i clusters", nrow(re_info)), "\n")
  # }
  # Filter
  if (!is.null(filt_function)) {
    features_included = filt_function(train_data, train_labels, filt_args)
    train_data = train_data[, features_included]
    test_data = test_data[, features_included]
    cat(sprintf("Filtered features down to %i", length(features_included)), "\n")
  }
  # Fit
  fit = cv_function(train_data, train_labels, cv_args)
  # Predict
  preds = data.frame(
    "preds" = as.numeric(predict(fit, test_data)), #, s = "lambda.min"
    "age" = test_labels,
    row.names = rownames(test_data)
  )
  res = list(
    "fit" = fit, 
    "preds" = preds,
    "seed" = seed)
  return(res)
}

predict_from_coefs = function(coefs, data, untransform_age=F) {
  intercept_ind = which(rownames(coefs) == "(Intercept)")
  feats = rownames(coefs)[-intercept_ind]
  preds = as.matrix(data[, feats]) %*% coefs[feats, "Coef"]
  preds = preds + coefs[intercept_ind, "Coef"]
  if (untransform_age) {
    preds = ageTinv(preds, 20)
  }
  return(preds)
}

predict_plot = function(preds, meta, color_var = NULL, untransform_age=NULL) {
  if (!is.null(untransform_age)) {
    preds$age = ageTinv(preds$age, untransform_age)
    preds$preds = ageTinv(preds$preds, untransform_age)
  } else if (!is.null(untransform_age) && untransform_age == T) {
    stop("Error: untransform age changed syntax")
  }
  tmp = drop_na(preds)
  anno = sprintf("RMSE = %f\nMAE = %f\nr = %f", 
                 rmse(tmp$age, tmp$preds),
                 mae(tmp$age,tmp$preds),
                 cor(tmp$age, tmp$preds))
  xrange = range(tmp$age)
  yrange = range(tmp$preds)
  if (!is.null(color_var)) {
    if (color_var == "fold") {
      p = ggplot(preds, aes(x=age, y=preds, color=fold))+
        geom_point(size=0.4)+
        labs(x="Age", y="Predicted age")+
        geom_abline(intercept=-5, slope=1, linetype="dashed")+
        geom_abline(intercept=0, slope=1)+
        geom_abline(intercept=+5, slope=1, linetype="dashed")+
        annotate("text", x = xrange[1], y=yrange[2], label=anno, size=3, vjust=1, hjust=0)+
        labs(color = "Fold")+
        lims(x=xrange, y=yrange)
    } else {
      preds = merge(preds, meta, by=0) %>%
        dplyr::rename(age = age.x)
      print(nrow(preds))
      p = ggplot(preds, aes(x=age, y=preds, color=.data[[color_var]]))+
        geom_point(size=0.4)+
        labs(x="Age", y="Predicted age")+
        geom_abline(intercept=-5, slope=1, linetype="dashed")+
        geom_abline(intercept=0, slope=1)+
        geom_abline(intercept=+5, slope=1, linetype="dashed")+
        annotate("text", x = xrange[1], y=yrange[2], label=anno, size=3, vjust=1, hjust=0)+
        labs(color = color_var)+
        lims(x=xrange, y=yrange)
    }
  } else {
    p = ggplot(preds, aes(x=age, y=preds))+
      geom_point(size=0.4)+
      labs(x="Age", y="Predicted age")+
      geom_abline(intercept=-5, slope=1, linetype="dashed")+
      geom_abline(intercept=0, slope=1)+
      geom_abline(intercept=+5, slope=1, linetype="dashed")+
      annotate("text", x = xrange[1], y=yrange[2], label=anno, size=3, vjust=1, hjust=0)+
      lims(x=xrange, y=yrange)
  }
  return(p)
}

##### CV functions #####

cv_glmnet = function(data, labels, args=list()) {
  alphas = args$alphas
  fold_ind = createFolds(labels, k = 10, list=F)
  best_fit = NULL
  best_mse = +Inf
  best_alpha = alphas[1]
  for (a in alphas) {
    fit = cv.glmnet(data, labels, foldid = fold_ind, alpha=a)
    if (fit$cvm[fit$index[1]] < best_mse) {
      best_fit = fit
      best_mse = fit$cvm[fit$index[1]]
      best_alpha = a
    }
  }
  # Because glmnet is annoying and returns the whole package, I refit manually
  best_lambda = fit$lambda.min
  final_fit = glmnet(data, labels, lambda=best_lambda, alpha=best_alpha)
  final_fit$alpha = best_alpha
  return(final_fit)
}

cv_glmnet_group = function(data, labels, groups, args=list()) {
  groups = as.integer(as.factor(groups))
  alphas = args$alphas
  best_fit = NULL
  best_mse = +Inf
  best_alpha = alphas[1]
  for (a in alphas) {
    fit = cv.glmnet(data, labels, foldid = groups, alpha=a) # Group is fold index
    if (fit$cvm[fit$index[1]] < best_mse) {
      best_fit = fit
      best_mse = fit$cvm[fit$index[1]]
      best_alpha = a
    }
  }
  # Because glmnet is annoying and returns the whole package, I refit manually
  best_lambda = fit$lambda.min
  final_fit = glmnet(data, labels, lambda=best_lambda, alpha=best_alpha)
  final_fit$alpha = best_alpha
  return(final_fit)
}

cv_glmnet_group_weighted = function(data, labels, groups, args=list()) {
  groups = as.integer(as.factor(groups))
  alphas = args$alphas
  best_fit = NULL
  best_mse = +Inf
  best_alpha = alphas[1]
  for (a in alphas) {
    fit = cv.glmnet(data, labels, foldid = groups, alpha=a, weights=args$weights) # Group is fold index
    if (fit$cvm[fit$index[1]] < best_mse) {
      best_fit = fit
      best_mse = fit$cvm[fit$index[1]]
      best_alpha = a
    }
  }
  # Because glmnet is annoying and returns the whole package, I refit manually
  best_lambda = fit$lambda.min
  final_fit = glmnet(data, labels, lambda=best_lambda, alpha=best_alpha)
  final_fit$alpha = best_alpha
  return(final_fit)
}

cv_glmnet_loo = function(data, labels, args=list()) {
  alphas = args$alphas
  fold_ind = createFolds(labels, k = length(labels), list=F)
  best_fit = NULL
  best_mse = +Inf
  best_alpha = alphas[1]
  for (a in alphas) {
    fit = cv.glmnet(data, labels, foldid = fold_ind, alpha=a)
    if (fit$cvm[fit$index[1]] < best_mse) {
      best_fit = fit
      best_mse = fit$cvm[fit$index[1]]
      best_alpha = a
    }
  }
  # Because glmnet is annoying and returns the whole package, I refit manually
  best_lambda = fit$lambda.min
  final_fit = glmnet(data, labels, lambda=best_lambda, alpha=best_alpha)
  final_fit$alpha = best_alpha
  return(final_fit)
}

cv_lm = function(data, labels, args=NULL) {
  fit = lm(labels~., data=data)
  return(fit)
}

cv_ffs = function(data, labels, args=NULL) {
  lower = lm(labels ~ 1, data=data)
  upper = lm(labels ~ ., data=data)
  fit = stepAIC(lower, data=data, direction = "forward", scope=list(upper=upper,lower=lower), 
                steps=args$steps, trace=F)
  return(fit)
}

##### Filtering #####

filt_by_age_cor = function(data, labels, args=list()) {
  cors = cor(data, labels)
  return(which(abs(cors) > args$cor_th))
}

filt_by_age_cor_crct = function(data, labels, args=list()) {
  data_crct = args$data_crct[rownames(data), colnames(data)]
  cors = cor(data_crct, labels)
  return(which(abs(cors) > args$cor_th))
}

filt_by_min_n = function(data, labels, args=list()) {
  return(which(args$re_info[colnames(data), "n"] >= args$minn))
}

##### Clustering #####

find_feature_clusters = function(data, re_info) {
  cors = cor(data)
  dists = as.dist(1-cors)
  clust = hclust(dists)
  cut_heights = seq(0, 2, 0.02)
  tmp = data.frame()
  for (h in cut_heights) {
    treecut = cutree(clust, h=h)
    rdm_cut = sample(treecut)
    homog = homogeneity(re_info$class, treecut)
    rdm_homog = homogeneity(re_info$class, rdm_cut)
    clust_classes = data.frame(
      clust = treecut,
      class = re_info$class) %>%
      group_by(clust) %>%
      summarize(n = n(), distinct_classes = n_distinct(class))
    random_clust = data.frame(
      clust = rdm_cut,
      class = re_info$class) %>%
      group_by(clust) %>%
      summarize(n = n(), distinct_classes = n_distinct(class))
    tmp = rbind(tmp, c(
      h, length(unique(treecut)), 
      homog, rdm_homog,
      mean(clust_classes$n), 
      mean(clust_classes$distinct_classes),
      mean(random_clust$distinct_classes)))
  }
  colnames(tmp) = c("height", "featn", "homogeneity", "rdm_homogeneity",  
                    "mean_re_per_clust",  "mean_classes_per_clust", "mean_classes_per_rdm_clust")
  out = list(
    tree = clust,
    cuts = tmp
  )
  return(out)
}

collapse_to_clusters = function(data, clusts, rows_are_features=F, na.rm=T) {
  uniq_clusts = sort(setdiff(unique(clusts), NA))
  data_clust = list()
  if (!rows_are_features) {
    for (c in uniq_clusts) {
      cols = which(!is.na(clusts) & clusts == c)
      data_clust[[c]] = rowMeans(data[cols], na.rm=na.rm)
    }
    data_clust = do.call(cbind, data_clust)
    rownames(data_clust) = rownames(data)
    colnames(data_clust) = uniq_clusts
  } else {
    for (c in uniq_clusts) {
      rows = which(!is.na(clusts) & clusts == c)
      data_clust[[c]] = colMeans(data[rows, ], na.rm=na.rm)
    }
    data_clust = do.call(rbind, data_clust)
    colnames(data_clust) = colnames(data)
    rownames(data_clust) = uniq_clusts
  }
  if (is.data.frame(data)) {
    data_clust = data.frame(data_clust)
  }
  return(data_clust)
}

cluster_treecut = function(data, re_info=NULL, args) {
  cors = cor(data)
  dists = as.dist(1-cors)
  tree = hclust(dists)
  clusts = cutree(tree, h=args$h)
  return(paste0("c", clusts))
}

cluster_min_n = function(data, re_info, args) {
  cors = cor(data)
  dists = as.dist(1-cors)
  tree = hclust(dists)
  clusts = 1:ncol(data)
  leaf_list = list()
  ns = re_info$n
  for (i in 1:nrow(tree$merge)) {
    lefti = tree$merge[i, 1]
    righti = tree$merge[i, 2]
    if (lefti < 0) {
      left_leaves = -lefti
    } else {
      left_leaves = leaf_list[[lefti]]
    }
    if (righti < 0) {
      right_leaves = -righti
    } else {
      right_leaves = leaf_list[[righti]]
    }
    leaf_list[[i]] = c(left_leaves, right_leaves)
    if (tree$height[i] < args$maxh & (sum(ns[left_leaves]) < args$minn | sum(ns[right_leaves]) < args$minn)) {
      clusts[leaf_list[[i]]] = min(clusts[left_leaves], clusts[right_leaves])
    }
  }
  return(paste0("c", clusts))
}

cluster_min_n_by_class = function(data, re_info, args) {
  for (cl in unique(re_info$class)) {
    featind = which(re_info$class == cl)
    if (length(featind) > 1) {
      re_info[featind, "clust"] = paste0(cl, "_", cluster_min_n(data[, featind], re_info, args=list(minn=args$minn, maxh=args$maxh)))
    } else {
      re_info[featind, "clust"] = paste0(cl, "_c1")
    }
  }
  return(re_info$clust)
}

##### Nanopore #####

collect_basecalling_stats = function(nanopath, stats) {
  # Estimated n of bases from sequencing run
  ont_report = dir(paste0(nanopath, "/reports"), pattern="json", full.names=T)
  ont_report = fromJSON(file = ont_report)
  stats$est_bases = 0
  # Sum over multiple acquisitions
  for (i in 1:length(ont_report$acquisitions)) {
    tmp = as.numeric(ont_report$acquisitions[[i]]$acquisition_run_info$yield_summary$estimated_selected_bases)
    stats$est_bases = stats$est_bases + tmp
  }
  # Reads passing basecalling
  files = dir(paste0(nanopath, "/basecalls/pass/cat_bams"), 
              pattern="lengths", full.names=T)
  fnames = str_extract(files, "/([[:alnum:]]+).lengths", group=1)
  stats$pass_n = c()
  stats$rlengths = list()
  for (i in 1:length(files)) {
    stats$rlengths[[fnames[i]]] = as.integer(readLines(files[i]))
    stats$pass_n[fnames[i]] = sum(stats$rlengths[[fnames[i]]])
  }
  # Passing with no or wrong bc
  stats$pass_wrong_bc_tot = 0
  files = dir(paste0(nanopath, "/basecalls/pass/cat_bams/hide"), 
              pattern="lengths", full.names=T)
  fnames = str_extract(files, "/([[:alnum:]]+).lengths", group=1)
  for (i in 1:length(files)) {
    if (fnames[i] == "unclassified") stats$pass_no_bc_tot = sum(as.integer(readLines(files[i])))
    else stats$pass_wrong_bc_tot = stats$pass_wrong_bc_tot + sum(as.integer(readLines(files[i])))
  }
  # Reads failing basecalling
  files = dir(paste0(nanopath, "/basecalls/fail/cat_bams"), 
              pattern="lengths", full.names=T)
  fnames = str_extract(files, "/([[:alnum:]]+).lengths", group=1)
  stats$fail_n = c()
  for (i in 1:length(files)) {
    stats$fail_n[fnames[i]] = sum(as.integer(readLines(files[i])))
  }
  # Passing with no or wrong bc
  stats$fail_wrong_bc_tot = 0
  files = dir(paste0(nanopath, "/basecalls/fail/cat_bams/hide"), 
              pattern="lengths", full.names=T)
  fnames = str_extract(files, "/([[:alnum:]]+).lengths", group=1)
  for (i in 1:length(files)) {
    if (fnames[i] == "unclassified") stats$fail_no_bc_tot = sum(as.integer(readLines(files[i])))
    else stats$fail_wrong_bc_tot = stats$fail_wrong_bc_tot + sum(as.integer(readLines(files[i])))
  }
  
  stats$pass_right_bc_tot = sum(stats$pass_n)
  stats$fail_right_bc_tot = sum(stats$fail_n)
  stats$pass_tot = stats$pass_right_bc_tot + stats$pass_wrong_bc_tot + stats$pass_no_bc_tot
  stats$fail_tot = stats$fail_right_bc_tot + stats$fail_wrong_bc_tot + stats$fail_no_bc_tot
  return(stats)
}

collect_mapping_stats = function(nanopath, stats) {
  files = dir(paste0(nanopath, "/mapped"), pattern="lengths", full.names=T)
  fnames = str_extract(files, "/([[:alnum:]]+).lengths", group=1)
  stats$mapped_n = c()
  for (i in 1:length(files)) {
    stats$mapped_n[fnames[i]] = sum(as.integer(readLines(files[i])))
  }
  stats$mapped_tot = sum(stats$mapped_n)
  return(stats)
}

read_nanopore = function(nanopath, anno_re, meta, stats) {
  # Make re ranges
  anno_ranges = GRanges(
    seqnames = anno_re$chr,
    names = anno_re$fragment_id,
    ranges = IRanges(
      start = anno_re$start,
      end = anno_re$end))
  # Lists to store one item for each file
  data = list()
  cnames = c("chr", "start", "end", "mod_type", "score", "strand", 
             "cov", "perc", "unmod", "mod", "filt", "nocall", "alt")
  files = dir(paste0(nanopath, "/basemods"), pattern="bed", full.names=T)
  fnames = str_extract(files, "/([[:alnum:]]+).cpg.acc.bed", group=1)
  stats$anycov = c()
  stats$quant = c()
  stats$noquant = c()
  stats$subdel = c()
  stats$anycovRE = c()
  stats$quantRE = c()
  stats$noquantRE = c()
  stats$subdelRE = c()
  for (i in 1:length(files)) {
    # Load modbases file
    tmp = fread(files[i], data.table=F)
    tmp = tmp[,-c(7:9)]
    colnames(tmp) = cnames
    # Note that subs and dels are a thing
    # Create a few utility variables
    tmp = tmp %>%
      mutate(quant = unmod + mod + alt) %>%
      mutate(noquant = filt + nocall) %>%
      mutate(subdel = cov - quant - noquant)
    # Save some stats
    stats$anycov[fnames[i]] = sum(tmp$cov)
    stats$quant[fnames[i]] = sum(tmp$quant)
    stats$noquant[fnames[i]] = sum(tmp$noquant)
    stats$subdel[fnames[i]] = sum(tmp$subdel)
    # Global
    global = tmp %>%
      dplyr::select(cov, quant, noquant, subdel, unmod, mod, alt) %>%
      summarize_all(sum) %>%
      mutate(featureID = "GenomeWide", reID = "GenomeWide",
             class = "GenomeWide", superf = "GenomeWide", 
             fam = "GenomeWide", name = "GenomeWide", .before=cov)
    # Collapse to repeat instances
    cpg_ranges = GRanges(
      seqnames = tmp$chr,
      ranges = IRanges(
        start = tmp$start,
        end = tmp$end))
    overlaps = findOverlaps(cpg_ranges, anno_ranges)
    overlaps = data.frame(overlaps)
    overlaps$fragment_id=anno_re$fragment_id[overlaps$subjectHits]
    tmp[overlaps$queryHits, "fragment_id"] = overlaps$fragment_id
    tmp = merge(anno_re, tmp, by="fragment_id") %>%
      dplyr::select(featureID, reID, class, superf, fam, name, 
                    cov, quant, noquant, subdel, unmod, mod, alt)
    tmp = tmp %>%
      group_by(featureID, reID, class, superf, fam, name) %>%
      summarize_all(sum) %>%
      dplyr::filter(cov > 0)
    # Save some stats
    stats$anycovRE[fnames[i]] = sum(tmp$cov)
    stats$quantRE[fnames[i]] = sum(tmp$quant)
    stats$noquantRE[fnames[i]] = sum(tmp$noquant)
    stats$subdelRE[fnames[i]] = sum(tmp$subdel)
    tmp = rbind(tmp, global)
    data[[fnames[i]]] = tmp
  }
  # Make homogenous lists
  data2 = list()
  data2$fivemc = data2$hmc = list(
    met = list(),
    cov = list()
  )
  data2$subdel = list(
    subdel = list(),
    cov = list()
  )
  for (bc in names(data)) {
    # 5mc met and cov
    tmp = subset(data[[bc]], quant > 0)
    tmp[bc] = tmp$mod
    data2$fivemc$met[[bc]] = tmp[c("featureID", "reID", "class", "superf", 'fam', "name", bc)]
    tmp[bc] = tmp$quant
    data2$fivemc$cov[[bc]] = tmp[c("featureID", "reID", "class", "superf", 'fam', "name", bc)]
    # 5hmc met and cov
    tmp = subset(data[[bc]], quant > 0)
    tmp[bc] = tmp$alt
    data2$hmc$met[[bc]] = tmp[c("featureID", "reID", "class", "superf", 'fam', "name", bc)]
    tmp[bc] = tmp$quant
    data2$hmc$cov[[bc]] = tmp[c("featureID", "reID", "class", "superf", 'fam', "name", bc)]
    # Any difference met and cov
    tmp = subset(data[[bc]], quant + subdel > 0)
    tmp[bc] = tmp$subdel
    data2$subdel$subdel[[bc]] = tmp[c("featureID", "reID", "class", "superf", 'fam', "name", bc)]
    tmp[bc] = tmp$quant + tmp$subdel
    data2$subdel$cov[[bc]] = tmp[c("featureID", "reID", "class", "superf", 'fam', "name", bc)]
  }
  # Concatenate lists to tables and rename 
  name_conv = c("featureID", "reID", "class", "superf", 'fam', "name", meta$ID)
  names(name_conv) = c("featureID", "reID", "class", "superf", 'fam', "name", meta$Barcode)
  for (n1 in names(data2)) {
    for (n2 in names(data2[[n1]])) {
      data2[[n1]][[n2]] = Reduce(function(x, y) merge(x, y, by=c("featureID", "reID", "class", "superf", 'fam', "name"), all=T), data2[[n1]][[n2]])
      colnames(data2[[n1]][[n2]]) = name_conv[colnames(data2[[n1]][[n2]])]
      data2[[n1]][[n2]][is.na(data2[[n1]][[n2]])] = 0
    }
  }
  for (n1 in names(data2)) {
    data2[[n1]]$re_info = data2[[n1]]$cov[, c("featureID", "reID", "class", "superf", 'fam', "name")]
    global_ind = which(data2[[n1]]$re_info$superf == "GenomeWide")
    data2[[n1]]$re_info = data2[[n1]]$re_info %>%
      bind_rows(dplyr::slice(., global_ind), .) %>%
      dplyr::filter(!row_number() %in% c(global_ind+1))
    rownames(data2[[n1]]$re_info) = c("global", data2[[n1]]$re_info$featureID[-1])
    for (n2 in setdiff(names(data2[[n1]]), "re_info")) {
      data2[[n1]][[n2]] = data2[[n1]][[n2]] %>%
        bind_rows(dplyr::slice(., global_ind), .) %>%
        dplyr::filter(!row_number() %in% c(global_ind+1))
      rownames(data2[[n1]][[n2]]) = rownames(data2[[n1]]$re_info)
      data2[[n1]][[n2]][c("featureID", "reID", "class", "superf", 'fam', "name")] = NULL
      data2[[n1]][[n2]] = data.frame(t(data2[[n1]][[n2]]))
    }
  }
  out = list(
    data = data2,
    stats = stats
  )
  return(out)
}

plot_losses = function(stats) {
  edges = data.frame(
    N1 = c("FlowCell", "FlowCell", "FlowCell",
           "PassBaseCall", "PassBaseCall", "PassBaseCall",
           "FailBaseCall", "FailBaseCall", "FailBaseCall",
           "PassRightBC", "PassRightBC",
           "Mapped",
           "AnyCov", "AnyCov", "AnyCov",
           "Quant", "Quant",
           "NoQuant", "NoQuant",
           "SubDel", "SubDel"),
    N2 = c("PassBaseCall", "FailBaseCall", "NotCalled",
           "PassRightBC", "PassWrongBC", "PassNoBC",
           "FailRightBC", "FailWrongBC", "FailNoBC",
           "Mapped", "UnmappedOrMMap",
           "AnyCov",
           "Quant", "SubDel", "NoQuant",
           "RE", "nonRE",
           "RE", "nonRE",
           "RE", "nonRE"),
    Value = c(stats$pass_tot, stats$fail_tot, stats$est_bases - stats$pass_tot - stats$fail_tot,
              stats$pass_right_bc_tot, stats$pass_wrong_bc_tot, stats$pass_no_bc_tot,
              stats$fail_right_bc_tot, stats$fail_wrong_bc_tot, stats$fail_no_bc_tot,
              stats$mapped_tot, stats$pass_right_bc_tot - stats$mapped_tot,
              100*sum(stats$anycov),
              100*sum(stats$quant), 100*sum(stats$subdel), 100*sum(stats$noquant),
              100*sum(stats$quantRE), 100*(sum(stats$quant)-sum(stats$quantRE)),
              100*sum(stats$noquantRE), 100*(sum(stats$noquant)-sum(stats$noquantRE)),
              100*sum(stats$subdelRE), 100*(sum(stats$subdel)-sum(stats$subdelRE)))
  )
  nodes = data.frame(
    ID = c("FlowCell", 
           "PassBaseCall", "FailBaseCall", "NotCalled",
           "PassRightBC", "PassWrongBC", "PassNoBC",
           "FailRightBC", "FailWrongBC", "FailNoBC",
           "Mapped", "UnmappedOrMMap",
           "AnyCov",
           "Quant", "SubDel", "NoQuant",
           "RE", "nonRE"),
    x = c(0, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 4, 5, 5, 5, 6, 6)
  )
  styles = list(
    "FlowCell" = list(col="yellow", srt=0, textpos=4),
    "PassBaseCall" = list(col="green", srt=0),
    "FailBaseCall" = list(col="red", srt=0),
    "NotCalled" = list(col="red", srt=0),
    "FailRightBC" = list(col="red", srt=0), 
    "FailWrongBC" = list(col="red", srt=0), 
    "FailNoBC" = list(col="red", srt=0),
    "PassRightBC" = list(col="green", srt=0), 
    "PassWrongBC" = list(col="red", srt=0), 
    "PassNoBC" = list(col="red", srt=0),
    "Mapped" = list(col="green", srt=0),
    "UnmappedOrMMap" = list(col="red", srt=0),
    "AnyCov" = list(col="green", srt=0),
    "Quant" = list(col="green", srt=0), 
    "SubDel" = list(col="yellow", srt=0), 
    "NoQuant" = list(col="red", srt=0),
    "RE" = list(col="green", srt=0, textpos=2),
    "nonRE" = list(col="red", srt=0, textpos=2)
  )
  rp = makeRiver(nodes = nodes, edges = edges, styles = styles)
  plot(rp)
}

#### OLDER ####

strat2 = function(met, cov, train_ind, theta=0) {
  # Partition
  train_met = met[train_ind,]
  train_cov = cov[train_ind,]
  train_labels = meta$age[train_ind]
  test_met = met[-train_ind,]
  test_cov = cov[-train_ind,]
  test_labels = meta$age[-train_ind]
  # Fit
  coefs = train_ML_clock(train_met, train_cov, train_labels)
  selected = coefs$padj < 0.001
  preds = loglike_pred(test_met[, selected], test_cov[, selected], coefs[selected, ], theta)
  res = data.frame(
    "ages" = test_labels,
    "preds" = as.vector(preds)
  )
  return(res)
}

##### ML CLOCK #####

# Log likelihood function for logit transformed
# i represents one sample
# coefi stores intercept, age coef, nuisance coefs
loglike_expit = function(age, meti, covi, coefs, nuisancei=NULL, theta=0) {
  w = 1 # Explore more later
  # w = log(covi)
  mm = t(replicate(nrow(coefs), c(1, age, unlist(nuisancei))))
  X = coefs * mm
  X = rowSums(X)
  L = meti*w*X - covi*w*log(1+exp(X))
  return(-sum(L))
}

loglike_expit_debug = function(age, meti, covi, coefs, nuisancei=NULL, theta=0) {
  w = 1 # Explore more later
  mm = t(replicate(nrow(coefs), c(1, age, unlist(nuisancei))))
  X = coefs * mm
  X = rowSums(X)
  L = meti*w*X - covi*w*log(1+exp(X))
  return(-L)
}

# Predict based on maximum likelihood
predict_ML_clock = function(met, cov, coefs, nuisance=NULL, theta=0, range=c(0, 120)) {
  if (!is(cov, "matrix") | !is(met, "matrix")) {
    stop("met and cov should be matrices")
  }
  preds = c()
  for (i in 1:nrow(met)) {
    opt = optimize(loglike_expit, interval = range, 
                   meti = met[i,rownames(coefs)], covi = cov[i,rownames(coefs)], 
                   coefs = coefs, nuisance=nuisance[i,], theta=theta)
    preds = c(preds, opt$minimum)
  }
  return(preds)
}

# Train coefficients for ML clock
train_ML_clock = function(met, cov, explanatory, ncoef=50, rank_by="r2") {
  met = as.matrix(met)
  cov = as.matrix(cov)
  frml = paste(colnames(explanatory), collapse="+")
  frml = as.formula(paste("Y~", frml))
  coefs = matrix(ncol=1+ncol(explanatory), nrow=ncol(cov))
  rownames(coefs) = colnames(cov)
  colnames(coefs) = c("intercept", colnames(explanatory))
  r2s = c()
  for (re in colnames(cov)) {
    Y = as.matrix(cbind(met[, re], cov[, re]-met[, re]))
    fit = glm(frml, data=explanatory, family = binomial())
    r2s[re] = with(summary(fit), 1 - deviance/null.deviance)
    coefs[re, ] = fit$coefficients
  }
  if (rank_by == "r2") {
    print(sprintf("Selecting %i features by %s", ncoef, rank_by))
    topN_coefs = order(-r2s)[1:ncoef]
  } else if(rank_by == "coef") {
    print(sprintf("Selecting %i features by %s", ncoef, rank_by))
    topN_coefs = order(-abs(coefs[,2]))[1:ncoef]
  }
  return(coefs[topN_coefs, ])
}

##### BIO AGE ESTIMATE #####

# Cost function for bio age optimization
cost = function(age, ms, coefs) {
  residuals = ms - rowSums(coefs[c("c0", "c1")] * matrix(c(1, age), ncol=2, nrow=nrow(coefs), byrow=T))
  sum(residuals^2)
}

# Cost function gradient for bio age optimization
costgr = function(age, ms, coefs) {
  residuals = ms - rowSums(coefs[c("c0", "c1")] * matrix(c(1, age), ncol=2, nrow=nrow(coefs), byrow=T))
  sum(-2*coefs["c1"]*residuals)
}

# Shift ages to estimate biological age
shift_ages = function(data, ages, coefs) {
  res = data.frame(
    cost_start = numeric(),
    cost_end = numeric(),
    iter = integer(),
    code = integer(),
    chrono_age = numeric(),
    bio_age = numeric()
  )
  for (i in 1:length(ages)) {
    res[i, "cost_start"] = cost(ages[i], ms=data[i, ], coefs = coefs)
    res[i, "chrono_age"] = ages[i]
    opt = optim(par=ages[i], fn=cost, gr = costgr, ms=data[i, ], coefs = coefs, 
                method="Brent", lower=-20, upper = 150)
    res[i, "cost_end"] = opt$value
    res[i, "iter"] = opt$counts[1]
    res[i, "code"] = opt$convergence
    res[i, "bio_age"] = opt$par
    print(i)
    # if (i == 100) break
  }
  rownames(res) = rownames(data)
  return(res)
}

collapse_to_clusters_old = function(data, re_info_ref, re_info_this=NULL) {
  if (is.null(re_info_this)) { # Then this re_info is the reference one
    re_info_this = re_info_ref
  } else {
    re_info_this = merge(re_info_this, re_info_ref[c("superf", "fam", "clust")], by=c("superf", "fam"))
  }
  data_clust = data.frame(t(data))
  data_clust$clust = re_info_this[rownames(data_clust), "clust"]
  data_clust = data_clust %>%
    group_by(clust) %>%
    summarize_all(mean) %>%
    column_to_rownames("clust")
  data_clust = t(data_clust)
  if (is.data.frame(data)) {
    data_clust = data.frame(data_clust)
  }
  out = list(
    data = data_clust,
    re_info = re_info_this
  )
  return(out)
}
