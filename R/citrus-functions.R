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
#' @import flowCore
#' @export
ReadFCSSet <- function (dataDirectory, fileList, fileSampleSize = 1000,
                        transformColumns = NULL, transformCofactor = 5,
                        scaleColumns = NULL, useChannelDescriptions = F,
                        readParameters = NULL, ...) {
  data = list()
  fileCounter = 1
  fileNames = c()
  fileChannelNames = list()
  fileReagentNames = list()
  addtlArgs = list(...)
  conditions = colnames(fileList)
  for (i in 1:length(conditions)) {
    cat(paste("Reading Condition ", conditions[i], "\n"))
    conditionData = list()
    fileChannelNames[[conditions[i]]] = list()
    fileReagentNames[[conditions[i]]] = list()
    for (fileName in fileList[, conditions[i]]) {
      fileNames[fileCounter] = fileName
      filePath = file.path(fileName)
      if (!file.exists(filePath)) {
        stop(paste("File", filePath, "not found."))
      }
      cat(paste("\tReading file ", fileName, "\n"))
      fcsFile = ReadFCS(filePath)
      fcsData = fcsFile@exprs
      parameterDescriptions = as.vector(pData(flowCore::parameters(fcsFile))$desc)
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
        cat(paste("\tSampling", fileSampleSize, "events.\n"))
        fcsData = fcsData[sort(sample(1:nrow(fcsData),
                                      fileSampleSize)), ]
      }
      conditionData[[fileName]] = fcsData
      print(dim(fcsData))
      if (useChannelDescriptions) {
        channelDescriptions = as.vector(pData(parameters(fcsFile))$desc)
        nchar(channelDescriptions) > 2
      }
    }
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
#' @param minimumClusterSizePercent Minimal size of cluster in proportion.
#' @param fileSampleSize The number of cells subsampled from each file.
#' @param nFolds Number of folds used for clustering
#' @param featureType A character vector of feature type (abundances, proportion).
#' @param seed A seed to reproduce the run.
#' @param pattern Regular expression used to select FCS files in the folder.
#' @return A list consisting of feature matrix, citrus clustering object and
#' combinedFCSSet object.
#' @import citrus
#' @export
runCITRUS = function(dataDir, selected.desc, fileList=NULL,
                            family = "classification",
                            minimumClusterSizePercent = 0.05,
                            fileSampleSize = 1000, nFolds = 1,
                            featureType = c("abundances"), seed=1,
                            pattern='fcs') {
  set.seed(seed)

  outputDir = file.path(dataDir, "citrusOutput")
  if (is.null(fileList)) {
    fileList = list.files(dataDir, pattern=pattern, full=T)
  }
  channels = GetParameters(fileList[1])
  colnames.map = hashmap::hashmap(channels$desc, channels$name)
  colnames.selected = sapply(selected.desc, function(x) colnames.map[[x]])
  clusteringColumns = colnames.selected
  transformColumns = colnames.selected
  transformCofactor = 5
  scaleColumns = colnames.selected
  fileList = data.frame(defaultCondition=fileList)

  combinedFCSSet = ReadFCSSet(dataDir, fileList, fileSampleSize,
                              transformColumns, transformCofactor)

  citrus.foldClustering = citrus.clusterAndMapFolds(combinedFCSSet,
                                                    clusteringColumns,
                                                    rep(NA, length(fileList)),
                                                    nFolds)

  conditions = colnames(fileList)[1]
  features = citrus.calculateFoldFeatureSet(citrus.foldClustering, combinedFCSSet,
                    featureType=featureType, conditions=conditions,
                    minimumClusterSizePercent=minimumClusterSizePercent)

  return(list(features=features, foldClustering=citrus.foldClustering,
              combinedFCS=combinedFCSSet))
}
