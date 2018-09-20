#' Run PhenoGraph on a FCS File
#'
#' Read FCS file and select columns to run PhenoGraph on.
#'
#' @param path Full path to the FCS file.
#' @param desc.selected An optional character vector to specify description of
#' column selected.
#' @param k Nearest-neighbor parameter of PhenoGraph.
#' @param seed An optional seed for reproducibility.
#' @return A Rphenograph object.
#' @export
RunPhenoGraph = function(path, desc.selected=NULL, k=30, TSNE=T, seed=1) {
  set.seed(1)
  fcs_data = flowCore::read.FCS(path)
  data = fcs_data@exprs
  colnames(data) = fcs_data@parameters@data$desc
  if(!is.null(desc.selected)) data = data[, desc.selected]
  phenograph_res = Rphenograph::Rphenograph(data=data, k=k)
  phenograph_tsne = NULL
  if(TSNE) {
    print('Running TSNE...')
    phenograph_tsne = PhenoGraphTSNE(data, desc.selected=desc.selected)
  }
  return(list(phenograph=phenograph_res, tsne=phenograph_tsne))
}

#' Run Rtsne on Data Matrix
#'
#' Run Rtsne on a data matrix on given columns.
#'
#' @param data A data matrix.
#' @param desc.selected Names of columns to run Rtsne on.
#' @return A Rtsne object.
PhenoGraphTSNE = function(data, desc.selected=NULL) {
  if(!is.null(desc.selected)) data = data[, desc.selected]
  tsne_out = Rtsne::Rtsne(as.matrix(data))
  return(tsne_out)
}

#' Plot TSNE using ggplot2
#'
#' Plot TSNE using ggplot2.
#' @param tsne_coord A data frame of first two dimension of tsne projection.
#' @param cluster A character/numeric vector cluster membership corresponding to
#' tsne_coord
#' @param sample.idx A numeric/logical vector indicating selected cells
#' @import ggplot2
#' @import dplyr
#' @export
PlotTSNE = function(tsne_coord, cluster, sample.idx=NULL, text=NULL) {
  tsne_coord = data.frame(tsne_coord)
  colnames(tsne_coord) = c('x', 'y')
  tsne_coord$cluster = cluster
  if(!is.null(sample.idx)) tsne_coord = tsne_coord[sample.idx,]
  tsne_cluster_median = tsne_coord %>%
    group_by(cluster) %>%
    summarize_all(median)

  p = ggplot(data=tsne_coord) +
    geom_point(aes(x=x, y=y, color=cluster), alpha=0.3) +
    geom_text(data=tsne_cluster_median, aes(x=x, y=y, label=cluster), size=7) +
    theme_classic() +
    theme(legend.position = 'none') +
    xlab("t-SNE dimension 1") +
    ylab("t-SNE dimension 2")
  if(!is.null(text)) {
    p = p + annotate('text', -Inf, Inf, label=text, hjust=0, vjust=1, size=10)
  }
  return(p)
}

#' Construct Abundance Matrix
#'
#' Construct abundance matrix given the cluster assignment from PhenoGraph
#'
#' @param data A data frame consisting of three columns: fileid, cluster and
#' outcome corresponding to each cell.
#' @export
ConstructAbundanceMatrix = function(data) {
  data %>%
    group_by(fileid) %>%
    mutate(total = n()) %>%
    group_by(fileid, outcome, cluster, total) %>%
    summarize(count = n()) %>%
    mutate(abundance = count/total) %>%
    select(fileid, outcome, cluster, abundance) %>%
    tidyr::spread(cluster, abundance, fill=0)
}

