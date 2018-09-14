#' Export combinedFCS as FCS files
#'
#' Export combinedFCS object from citrus into a FCS file with selected channels.
#'
#' Exported files will be named as sampled.fcs in the folder specified. The
#' file that each cells in the exported FCS are indexed and saved as fileID.csv
#' and corresponding file names of each index are saved as filenames.csv
#' @param clustering A list containing combinedFCS object output from citrus.
#' @param out.folder A character scalar of full path to output folder.
#' @param name A character vector containing selected channel name.
#' @param desc A character vector containing corresponding channel description.
#' @return None
#' @export
ExportSubsetFCS = function(clustering, out.folder, name, desc) {
  datafiles = clustering$combinedFCS$data[,name]
  fileID = clustering$combinedFCS$data[,'fileId']
  filenames = clustering$combinedFCS$fileNames

  write.table(fileID, file=paste0(out.folder, '/fileID.csv'), row.names=F)
  write.table(filenames, file=paste0(out.folder, '/filenames.csv'), row.names=F)

  meta = data.frame(name=name, desc=desc)
  meta$range <- apply(apply(datafiles,2,range),2,diff)
  meta$minRange <- apply(datafiles,2,min)
  meta$maxRange <- apply(datafiles,2,max)
  parameter = Biobase::AnnotatedDataFrame(meta)
  ff.concat = flowFrame(exprs = datafiles, parameter=parameter)

  write.FCS(ff.concat, file=paste0(out.folder, '/sampled.fcs'))
}

#' Get FCS Channel Name and Description
#'
#' Get the FCS channel name and corresponding description.
#' @param path A character vector of FCS file full paths.
#' @return A data frame that contains the name and description of the channels
#' @export
GetParameters = function(path) {
  fs = flowCore::read.flowSet(path)
  desc = fs[[1]]@parameters$desc
  name = fs[[1]]@parameters$name
  ret = data.frame(i=seq(length(desc)), name=name, desc=desc)
  return(ret)
}
