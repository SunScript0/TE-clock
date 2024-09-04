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

##### Clustering #####

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