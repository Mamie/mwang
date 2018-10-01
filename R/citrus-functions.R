#' Read an FCS file
#'
#' Read an FCS file while suppressing warnings.
#' @param filePath The full path to the FCS files.
#' @param ... Additional arguments to flowCore::read.FCS.
#' @return A flowCore::flowFrame object.
#' @export
ReadFCS <- function (filePath, ...) {
  addtlArgs = list(...)
  dataset = 1
  if ("dataset" %in% names(addtlArgs))
    dataset = addtlArgs[["dataset"]]
  which.lines = NULL
  if ("which.lines" %in% names(addtlArgs))
    which.lines = addtlArgs[["which.lines"]]
  fcs = tryCatch({
    suppressWarnings(flowCore::read.FCS(filePath, dataset = dataset,
                              which.lines = which.lines,
                              truncate_max_range = F))
  }, error = function(e) {
    if (grepl("Please set argument 'emptyValue' as", e$message)) {
      suppressWarnings(flowCore::read.FCS(filePath, dataset = dataset,
                                which.lines = which.lines, emptyValue = F,
                                truncate_max_range = F))
    }
    else {
      stop(e$message)
    }
  })
  return(fcs)
}

#' @inherit citrus::citrus.readFCSSet
#' @export
ReadFCSSet <- function (dataDirectory, fileList, fileSampleSize = 1000,
                        transformColumns = NULL, transformCofactor = 5,
                        scaleColumns = NULL, useChannelDescriptions = F,
                        readParameters = NULL, verbose = TRUE, ...) {
  data = list()
  fileCounter = 1
  fileNames = c()
  fileChannelNames = list()
  fileReagentNames = list()
  addtlArgs = list(...)
  conditions = colnames(fileList)
  for (i in 1:length(conditions)) {
    if (verbose) cat(paste("Reading Condition ", conditions[i], "\n"))
    conditionData = list()
    fileChannelNames[[conditions[i]]] = list()
    fileReagentNames[[conditions[i]]] = list()
    p <- dplyr::progress_estimated(length(fileList[, conditions[i]]))
    cat("Reading in FCS files.\n")
    for (fileName in fileList[, conditions[i]]) {
      fileNames[fileCounter] = fileName
      filePath = file.path(fileName)
      if (!file.exists(filePath)) {
        stop(paste("File", filePath, "not found."))
      }
      if(verbose) cat(paste("\tReading file ", fileName, "\n"))
      fcsFile = ReadFCS(filePath)
      fcsData = fcsFile@exprs
      parameterDescriptions = as.vector(flowCore::pData(flowCore::parameters(fcsFile))$desc)
      parameterNames = flowCore::colnames(fcsFile)
      invalidDescriptions = unname(which(sapply(parameterDescriptions,
                                                nchar) < 3 | is.na(parameterDescriptions)))
      parameterDescriptions[invalidDescriptions] = parameterNames[invalidDescriptions]
      if (useChannelDescriptions) {
        colnames(fcsData) = parameterDescriptions
      }
      if (!is.null(readParameters)) {
        fcsData = fcsData[, readParameters]
      }
      fileChannelNames[[conditions[i]]][[fileName]] = parameterNames
      fileReagentNames[[conditions[i]]][[fileName]] = parameterDescriptions
      fcsData = cbind(fcsData, fileEventNumber = 1:nrow(fcsData),
                      fileId = fileCounter)
      fileCounter = fileCounter + 1
      if ((!is.null(fileSampleSize)) && (fileSampleSize <
                                         nrow(fcsData))) {
        if (verbose) cat(paste("\tSampling", fileSampleSize, "events.\n"))
        fcsData = fcsData[sort(sample(1:nrow(fcsData),
                                      fileSampleSize)), ]
      }
      conditionData[[fileName]] = fcsData
      if (verbose) cat(paste("\t FCS file size is", str(dim(fcsData)), '.\n'))
      if (useChannelDescriptions) {
        channelDescriptions = as.vector(flowCore::pData(parameters(fcsFile))$desc)
        nchar(channelDescriptions) > 2
      }
      cat(capture.output(p$tick()))
    }
    cat('\n')
    data[[conditions[i]]] = do.call("rbind", conditionData)
    rm(conditionData)
    gc()
  }
  data = do.call("rbind", data)
  results = list(fileIds = matrix(1:(fileCounter - 1), ncol = length(conditions),
                                  dimnames = list(c(), conditions)),
                 fileNames = fileNames, fileChannelNames = fileChannelNames,
                 fileReagentNames = fileReagentNames)
  if (!is.null(transformColumns)) {
    if (any(!is.numeric(transformColumns))) {
      containedCols = setdiff(transformColumns, colnames(data))
      if (length(containedCols) > 0) {
        stop(paste("Transform cols", paste(containedCols, collapse = ", "),
                   "not found. Valid channel names:",
                   paste(colnames(data), collapse = ", ")))
      }
    }
    data[, transformColumns] = asinh(data[, transformColumns]/transformCofactor)
    results$transformColumns = transformColumns
    results$transformCofactor = transformCofactor
  }
  if (!is.null(scaleColumns)) {
    if (any(!is.numeric(scaleColumns))) {
      containedCols = setdiff(scaleColumns, colnames(data))
      if (length(containedCols) > 0) {
        stop(paste("Scale cols", paste(containedCols, collapse = ", "),
                   "not found. Valid channel names:",
                   paste(colnames(data), collapse = ", ")))
      }
    }
    results$scaleColumns = scaleColumns
    results$scaleColumns.mean = apply(data[, scaleColumns], 2, mean)
    results$scaleColumns.SD = apply(data[, scaleColumns], 2, sd)
    data[, scaleColumns] = apply(data[, scaleColumns], 2, scale)
  }
  results$data = data
  class(results) = "citrus.combinedFCSSet"
  return(results)
}

#' Run CITRUS
#'
#' Run CITRUS on selected channels
#' @param dataDir Full path to the folder containing FCS files.
#' @param selected.desc A character vector of selected channel description.
#' @param fileList A list of full path to FCS files.
#' @param family A character scalar of model type.
#' @param transformCofactor An integer for arcsinh transform channel; default 5
#' @param minimumClusterSizePercent Minimal size of cluster in proportion.
#' @param fileSampleSize The number of cells subsampled from each file.
#' @param nFolds Number of folds used for clustering
#' @param featureType A character vector of feature type (abundances, proportion).
#' @param seed A seed to reproduce the run.
#' @param pattern Regular expression used to select FCS files in the folder.
#' @param labels Optional vector of label used for balancing clustering folds
#' @return A list consisting of feature matrix, citrus clustering object and
#' combinedFCSSet object.
#' @import citrus
#' @export
runCITRUS = function(dataDir, selected.desc, fileList=NULL,
                            family = "classification", transformCofactor = 5,
                            minClusterSizePercent = 0.05,
                            fileSampleSize = 1000, nFolds = 1,
                            featureType = c("abundances"), seed=1,
                            pattern='fcs', labels=NULL) {
  set.seed(seed)

  outputDir = file.path(dataDir, "citrusOutput")
  if (is.null(fileList)) {
    fileList = list.files(dataDir, pattern=pattern, full=T)
  }
  channels = GetParameters(fileList[1])
  colnames.map = hashmap::hashmap(channels$desc, channels$name)
  colnames.selected = sapply(selected.desc, function(x) colnames.map[[x]])
  fileList = data.frame(defaultCondition=fileList)

  combinedFCSSet <- ReadFCSSet(dataDir, fileList, fileSampleSize = fileSampleSize,
                              transformColumns = colnames.selected,
                              transformCofactor = transformCofactor,
                              scaleColumns = colnames.selected, verbose=FALSE)

  citrus.foldClustering <- citrus.clusterAndMapFolds(combinedFCSSet,
                                                    clusteringColumns = colnames.selected,
                                                    labels = labels,
                                                    nFolds = nFolds)
  conditions = colnames(fileList)[1]
  features = citrus.calculateFoldFeatureSet(citrus.foldClustering,
                                            combinedFCSSet,
                                            featureType=featureType,
                                            conditions=conditions,
                                            minimumClusterSizePercent=minClusterSizePercent)

  return(list(features=features,
              foldClustering=citrus.foldClustering,
              combinedFCS=combinedFCSSet))
}

#' CITRUS regression
#'
#' Run CITRUS regression function.
#'
#' @param clustering A list object returned by RunCITRUS.
#' @param outcome A vector of response.
#' @return A regression result object.
#' @export
CITRUSRegression <- function(clustering, outcome, seed = 1) {
  set.seed(seed)
  suppressWarnings({
    citrus.regressionResults <- mclapply(c("glmnet"),
                                        citrus.endpointRegress,
                                        citrus.foldFeatureSet=clustering$features,
                                        labels=outcome,
                                        family="classification")
  })
  return(citrus.regressionResults)
}

#' Plot CITRUS hierarchy
#'
#' Plot the hierarchy used in CITRUS colored by stratefying clusters.
#'
#' @param clustering A list object returned by RunCITRUS.
#' @param regressionRes An object returned by CITRUSRegression.
#' @param seed An optional seed for reproducibility.
#' @return A ggplot object of the hierarchy.
#' @export
PlotHierarchy <- function(clustering, regressionRes, seed=1) {
  set.seed(seed)

  diff_cluster <- data.frame(cluster=integer(), selected=logical())
  cluster_min <- regressionRes[[1]]$differentialFeatures$cv.min$clusters
  if (length(cluster_min)) {
    diff_cluster <- rbind(data.frame(cluster = cluster_min, selected = T))
  }
  graph.list <- citrus::citrus.createHierarchyGraph(clustering$foldClustering$allClustering,
                                                   clustering$features$allLargeEnoughClusters)
  g <- graph.list$graph

  cluster.attributes <- data.frame(cluster=as.numeric(as.character(vertex_attr(g, name='label')))) %>%
    dplyr::left_join(diff_cluster, by = "cluster") %>%
    tidyr::replace_na(list(selected=F))
  assignment <- clustering$foldClustering$allClustering$clusterMembership
  n <- c()
  for (i in seq(nrow(cluster.attributes))) {
    cluster.data <- clustering$combinedFCS$data[assignment[[cluster.attributes$cluster[i]]],]
    n <- c(n, nrow(cluster.data))
  }
  cluster.attributes$size <- n
  g <- set.vertex.attribute(g, "alpha", value=sapply(cluster.attributes$selected,
                                                    function(x) ifelse(x, 1, 0.3)))
  g <- set.vertex.attribute(g, "size", value=log2(cluster.attributes$size))
  net <- intergraph::asNetwork(g)
  p <- GGally::ggnet2(net, label='label', alpha="alpha", size='size',
                      edge.size=0.1, label.size=2.5, label.color='#616568',
                      palette = "Set2") +
    guides(size=FALSE, color=guide_legend(title=NULL)) +
    theme(legend.position='none')
  return(list(p=p, cluster.attributes = cluster.attributes))
}

#' Plot CITRUS Cluster MFI Heatmap
#'
#' Plot heatmap of mean fluroresence intensity of CITRUS clusters.
#' @param clustering A list object returned by RunCITRUS.
#' @param cluster.attributes A data frame returned by PlotHierarchy.
#' @param name Name of selected channels.
#' @param desc Description of selected channels.
#' @return A pheatmap object.
#' @export
PlotClusterMFI = function(clustering, cluster.attributes, name, desc, ...) {
  cluster.mean <- c()
  assignment = clustering$foldClustering$allClustering$clusterMembership
  clusters = cluster.attributes$cluster %>% as.character %>% as.numeric
  for (i in seq(clusters)) {
    cluster.data = clustering$combinedFCS$data[assignment[[clusters[i]]],]
    cluster.mean <- rbind(cluster.mean,
                          unlist(apply(cluster.data[,name], 2, mean)))
  }
  colnames(cluster.mean) = desc
  rownames(cluster.mean) = unlist(clusters)
  head(cluster.mean)

  cluster.attributes %<>% dplyr::mutate(label=ifelse(selected, '*', '')) %>%
    tidyr::unite(clusterlabeled, c('cluster', 'label'), sep='', remove=F)
  p.data <- data.matrix(cluster.mean)

  #for(i in seq(10)) {
  #  p.data = t(apply(p.data, 1, scale))
  #  p.data = apply(p.data, 2, scale)
  #}
  rownames(p.data) = rownames(cluster.mean)
  colnames(p.data) = colnames(cluster.mean)
  annotation = data.frame(selected=cluster.attributes$selected %>%
                            as.character)
  rownames(annotation) = cluster.attributes$cluster
  options(repr.plot.height=5.5, repr.plot.width=5)
  p = pheatmap::pheatmap(p.data, annotation_row=annotation, annotation_colors = list(
    selected = c('FALSE' = "white", 'TRUE' = "firebrick")), ...)
  return(p)
}

#' Plot Results of CITRUS regression
#'
#' Plot results of CITRUS regression including effect size of the
#' difference in abundance between groups, boxplot of abundance
#' between groups and distribution of differential clusters.
#' @param clustering A list object returned by RunCITRUS.
#' @param cluster.attributes A object returned by PlotHierarchy.
#' @param outcome A vector of binary outcome.
#' @param name Name of selected channels.
#' @param desc Description of selected channels.
#' @return A ggplot object of effect size.
#' @export
PlotRes = function(clustering, cluster.attributes, outcomes, name, desc, ylab='ES') {
  abundance = data.frame(clustering$features$allFeatures)
  colnames(abundance) = sapply(colnames(abundance),
                               function(x) strsplit(x, '\\.')[[1]][2])
  abundance$ID = seq(nrow(abundance))
  abundance$outcome = outcomes
  diff_clusters =  cluster.attributes[cluster.attributes$selected,]$cluster
  abundance %<>%
    tidyr::gather(cluster, level, -c(ID, outcome)) %>%
    dplyr::mutate(cluster = as.numeric(cluster))
  p.effsize = abundance %>%
    dplyr::filter(cluster %in% diff_clusters) %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(level=as.numeric(level)) %>%
    dplyr::summarize(d = -effsize::cohen.d(level, outcome)$estimate,
                     lwr = -effsize::cohen.d(level, outcome)$conf.int[1],
                     upr = -effsize::cohen.d(level, outcome)$conf.int[2]) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(cluster.attributes)
  ordering.clusters = p.effsize$cluster[order(p.effsize$d)]
  p.effsize = p.effsize %>%
    dplyr::mutate(cluster=factor(cluster, levels=ordering.clusters)) %>%
    na.omit() %>%
    ggplot(data=.) +
    geom_point(aes(x=cluster, y=d)) +
    geom_segment(aes(x=cluster, xend=cluster,
                     y=lwr, yend=upr), size=0.2) +
    coord_flip() +
    geom_hline(aes(yintercept=0), linetype='dashed', color='gray', size=0.5) +
    theme_classic() + ylab(ylab)

  p.boxplot = abundance %>%
    dplyr::filter(cluster %in% diff_clusters) %>%
    dplyr::mutate(cluster = factor(cluster, levels=rev(ordering.clusters))) %>%
    dplyr::mutate(level = as.numeric(level)) %>%
    ggplot(data=.) +
    geom_boxplot(aes(x=outcome, y=level), width=0.3) +
    geom_jitter(aes(x=outcome, y=level), color='steelblue', size=0.3)+
    ylab('abundance') +
    facet_wrap(~cluster, ncol=1) +
    theme_classic() +
    theme(strip.background=element_blank(), axis.title.y=element_blank()) +
    coord_flip()

  p.dist = PlotDist(clustering, name, desc, rev(ordering.clusters))
  return(list(p.effsize=p.effsize, p.boxplot=p.boxplot, p.dist=p.dist))
}

#' Plot Marker Distribution
#'
#' Plot marker distribution of each cluster.
#'
#' @param clustering A list object from RunCITRUS
#' @param name Name of selected channels.
#' @param desc Description of selected channels.
#' @param mincluster Clusters
PlotDist = function(clustering, name, desc, diff_cluster) {
  marker.num = length(name)
  signif.feat = as.numeric(as.character(diff_cluster))
  assignment = clustering$foldClustering$allClustering$clusterMembership

  background = data.frame(clustering$combinedFCS$data[,name])
  colnames(background) =  desc
  background = background %>%
    tidyr::gather(marker, level)
  p = list()
  for (i in seq(length(diff_cluster))) {
    cluster.df = data.frame(clustering$combinedFCS$data[assignment[[diff_cluster[i]]], name])
    colnames(cluster.df) = desc
    p[[i]] = cluster.df %>%
      tidyr::gather(marker, level) %>%
      ggplot(data=.) +
      geom_density(data=background, aes(x=level), fill='gray', alpha=0.5, size=0.) +
      geom_density(aes(x=level), fill='blue', alpha=0.5, size=0.1) +
      facet_wrap(~marker, ncol=marker.num, scale='free_x') +
      theme_classic() +
      theme(strip.background=element_blank(),
            axis.title.x=element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank()) +
      ylab(diff_cluster[i])
  }
  p = cowplot::plot_grid(plotlist=p, ncol=1)
  return(p)
}
