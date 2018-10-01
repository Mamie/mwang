## ----functions, include=FALSE--------------------------------------------
library(mwang)
filename_to_pid <- function(x) {
  x <- stringr::str_match(x, 'UPCC[0-9_-]+')[1]
  if (is.na(x)) return(NA)
  x <- strsplit(x, split='-|_', fixed=F)[[1]]
  trial <- x[1]
  pid <- ifelse(nchar(x[2]) == 2, x[2], as.character(as.numeric(x[3])))
  return(paste(trial, pid, sep='-'))
}

## ----file-path, message=FALSE, warning=FALSE-----------------------------
outcome_path <- "/Volumes/pdcs/szmamie/data/20180929_CLL_reanalysis/CLL_outcomes.csv"
TCD_dir_path <- "/Volumes/pdcs/szmamie/data/metaAnalysis/Apheresis/CD8/TCD/"
TCD_file_list <- list.files(file.path(TCD_dir_path, "8_20160823_ALL_CLL/"), pattern="fcs", full=T)
TCD_file_list_pid <- unname(sapply(TCD_file_list, filename_to_pid))
outcome_data <- readr::read_csv(outcome_path)
patient_selected <- outcome_data$lulian_outcome %in% c('CR', 'NR')
TCD_file_selected <- TCD_file_list[TCD_file_list_pid %in% outcome_data$pid[patient_selected]]
kableExtra::kable_styling(kableExtra::kable(outcome_data[patient_selected,]))

## ----load-data, message=FALSE, warning=FALSE-----------------------------
TCD_data_example <- mwang::ReadFCSSet(TCD_dir_path, data.frame(default=TCD_file_selected[1]))
channel_names <- TCD_data_example$fileChannelNames$default[[1]]
channel_desc <- TCD_data_example$fileReagentNames$default[[1]]
kableExtra::kable_styling(kableExtra::kable(cbind(seq(length(channel_desc)), channel_names, channel_desc)))
# select channels of interested by index
channel_selected_idx <- c(4, 5, 6, 7, 10, 12, 13, 15)
channel_desc_selected <- channel_desc[channel_selected_idx]
channel_names_selected <- channel_names[channel_selected_idx]

## ----run-citrus, message=FALSE, warning=FALSE----------------------------
clustering <- mwang::runCITRUS(dataDir = TCD_dir_path, 
                               selected.desc = channel_desc_selected, 
                               fileList = TCD_file_selected
                               )

## ----citrusRegression, message=FALSE, warning=FALSE----------------------
outcomes <- outcome_data$lulian_outcome[patient_selected]
model <- mwang::CITRUSRegression(clustering, outcomes)

## ---- fig.width=4, fig.height=4------------------------------------------
hierarchy <- mwang::PlotHierarchy(clustering, model)
hierarchy$p

## ----citrus-cluster, fig.width=12, fig.height=3--------------------------
res_plotlist = mwang::PlotRes(clustering, hierarchy$cluster.attributes, outcomes, channel_names_selected, 
        channel_desc_selected, ylab='ES')
cowplot::plot_grid(plotlist=res_plotlist, nrow=1, rel_widths=c(1, 1, 6)) 

