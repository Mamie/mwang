#' Run Lasso Regression
#'
#' Run Lasso Regression on abundance matrix provided.
#' @param X A data matrix consists of the predictors.
#' @param y A numeric/factor vector of response.
#' @return A list of cross validation and model object
#' @export
RunLasso = function(X, y, seed=1, ...) {
  set.seed(seed)
  suppressWarnings({
    cv.out = glmnet::cv.glmnet(X, y, family='binomial', nfolds=nrow(X), alpha=1, ...)
    })
  mod = glmnet::glmnet(X, y, family='binomial', alpha=1, lambda=cv.out$lambda.1min, ...)
  res = list(cv=cv.out, mod=mod)
  class(res) = c(class(res), 'RunLasso')
  return(res)
}

#' Plot Lasso Coefficient Trace
#'
#' Plot Lasso coefficient traces.
#' @param res A RunLasso object.
#' @export
plot.RunLasso = function(res) {
  par(mfrow=c(1, 2))
  plot(res$cv)
  plot(res$mod, "lambda", col = RColorBrewer::brewer.pal(12, "Set3"), lwd = sqrt(3))
}


#' Plot Marker Density
#'
#' Plot marker density plot
#' @param df A data frame with columns marker and corresponding cluster.
#' @param coefs The coefficient object from glmnet.
#' @return A ggplot object of the marker distribution for clusters with non-zero
#' coefficients
#' @export
PlotMarkerDensity = function(df, coefs) {
  background = df %>%
    select(-c(cluster)) %>%
    tidyr::gather(marker, level)

  marker.num = length(unique(background$marker))
  clusters = df %>%
    tidyr::gather(marker, level, -cluster)
  signif_feat = rownames(coefs)[coefs[,1]!=0][-1]
  p = list()
  for (i in seq(length(signif_feat))) {
    p[[i]] = clusters %>%
      dplyr::filter(cluster == signif_feat[i]) %>%
      ggplot(data=.) +
      geom_density(data=background, aes(x=level),
                   fill='gray', alpha=0.5, size=0.) +
      geom_density(aes(x=level), fill='blue', alpha=0.5, size=0.1) +
      facet_wrap(~marker, ncol=marker.num, scale='free_x') +
      theme_classic() + theme(strip.background=element_blank(),
                              axis.title.x = element_blank(),
                              axis.text = element_blank(),
                              axis.ticks = element_blank()) +
      ylab(signif_feat[i])
  }
  p.marker = cowplot::plot_grid(plotlist=p, ncol=1)
  return(p.marker)
}

#' Plot Abundance Boxplot
#'
#' Plot abundance boxplot between outcomes.
#'
#' @param abundance A data matrix that contains the abundance of each cluster and
#' corresponding outcome.
#' @return A ggplot object
#' @export
PlotAbundaceBoxplot = function(abundance, coefs) {
  signif_feat = rownames(coefs)[coefs[,1]!=0][-1]
  p.boxplot = abundance %>% ungroup %>% select(c('outcome', signif_feat)) %>%
    tidyr::gather(feature, level, -outcome) %>%
    mutate(feature = as.factor(as.numeric(feature))) %>%
    ggplot(data=.) +
    geom_boxplot(aes(outcome, level), width=0.4, outlier.size=0.8) +
    geom_jitter(aes(outcome, level), size=0.3, color='steelblue') +
    facet_wrap(~feature, ncol=1, scale='free_x') +
    coord_flip() + ylab('abundance') +
    theme_classic() +
    theme(strip.background=element_blank(),
          axis.title.y=element_blank())
  return(p.boxplot)
}
