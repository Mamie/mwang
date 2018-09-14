#' Extract Excel Sheet Names
#'
#' Extract Excel sheet names allowing regular expression for matching.
#' @param path A character scalar of full path name of the excel sheet.
#' @param pattern An optional regular expression. Only sheet names matching the
#' pattern will be returned. If not specified, all sheets are returned.
#' @return A character vector of sheet names matching the regular expression
#' @examples
#' \dontrun{FindSheet('path/to/xlsx', '^a-l')}
#' @export
ExtractSheetNames = function(path, pattern='') {
  sheet.names = readxl::excel_sheets(path)
  sheet.names[grepl(pattern, sheet.names)]
}

#' Read Excel Sheets into a List
#'
#' Read specified sheets of Excel files into a list.
#' @param path A character scalar of full path name of the excel sheet.
#' @param sheets A character vector of sheet names.
#' @param colnames A logical scalar to indicate if first row contain column names
#' @return A list of data frame of corresponding sheets
#' @export
ReadExcelSheet = function(path, sheets, colnames=F) {
  purrr::map(sheets, function(x) {
    raw = readxl::read_excel(path, sheet=x, trim_ws=T, col_names=colnames)
    stats::na.omit(raw)
  })
}

#' Process qPCR Data Frame List
#'
#' Process qPCR data frame in the list into the same format.
#' @param data A list of data frame containing qPCR measurements
#' @param qPCR.idx A numeric vector of size indicating the index of raw
#' and corrected qPCR values
#' @param select.idx A numeric vector indicating the indices of columns to be
#' selected other than those in `qPCR.idx`
#' @return A list of data frame with corrected qPCR values and selected columns
#' @export
ProcessqPCRList = function(data, qPCR.idx, select.idx) {
  header = data[[1]][1,]
  for (i in seq(length(data))) {
    raw = qPCR.idx[1]
    correct = qPCR.idx[2]
    data[[i]] = data[[i]][-1,]

    if(dim(data[[i]])[2] >= max(qPCR.idx)) {
      colnames(data[[i]]) = header
      for (j in seq(dim(data[[i]])[1])) {
        if (!grepl('lot', tolower(data[[i]][j, correct]))) {
          data[[i]][j, raw] = data[[i]][j, correct]
        }
      }
    } else {
      colnames(data[[i]]) = header[seq(dim(data[[i]])[2])]
    }
    data[[i]] = data[[i]][,c(select.idx, raw)]
  }
  data
}
