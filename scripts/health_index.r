#' Import the table of "good" and "bad" bacteria from Abass et al. 2021 metaanalysis paper.
#' amnona/paper-metaanalysis.
#' @param downlistpath path for tsv file of "bad" bacteria
#' @param uplistpath path for tsv file of "good" bacteria
#' @return data.frame with Feature_ID and dir (for bad and good bacteria)
#' @export import_nsf
#' @author Rotem Hadar
import_nsf <- function(downlistpath, uplistpath){
  location_prefix <- "https://raw.githubusercontent.com/amnona/paper-metaanalysis/main/ratios/nonspecific/"
  if(missing(downlistpath)) downlistpath <- paste0(location_prefix,"nonspecific-down_feature.txt")
  if(missing(uplistpath)) uplistpath <- paste0(location_prefix,"nonspecific-up_feature.txt")
  nsd <- readr::read_tsv(downlistpath)
  nsd <- dplyr::mutate(nsd, dir = "down")
  nsu <- readr::read_tsv(uplistpath)
  nsu <- dplyr::mutate(nsu, dir = "up")
  ns <- rbind(nsd, nsu)
  return(ns)
}

#' Calculate the Binary health index for each sample in an experiment.
#' @param exp data.frame of ASV experiment. ASV are rows and samples are the columns.
#' @param ns data.frame of "good" and "bad" bacteria with "dir".
#' @param thresh float. Threshold to consider level of bacteria absent/not.
#' @return data.frame. Binary table with score for each sample.
#' @export binary_health_index
binary_health_index <- function(exp, ns, thresh = 0){
  if(missing(ns)) ns <- import_nsf()
  # filter only relevant ASV, and pivot_longer
  exp <- 
    exp %>% data.frame() %>% 
    tibble::rownames_to_column("Feature_ID") %>% 
    dplyr::right_join(ns) %>% 
    tidyr::pivot_longer(contains("DB"), values_to = "reads", names_to = "sampleid") 
  # calculate HI
  res <- exp %>% 
    dplyr::mutate(exist = reads > thresh) %>% 
    dplyr::group_by(sampleid, dir) %>% dplyr::summarise(number = sum(exist, na.rm = T)) %>% 
    tidyr::pivot_wider(id_cols = sampleid, names_from = dir, values_from = number) %>% 
    dplyr::transmute(binary_HI = log2(down+.01)/log2(up+.01))
  res
}


#' Calculate frequency based Health Index for each sample in an ASV table.
#' @param exp data.frame of ASV exeriment. ASV are rows and samples are the columns.
#' @param ns data.frame of "good" and "bad" bacteria with "dir".
#' @param thresh floa. Threshold to consider level of bacteria absent/not.
#' @return data.frame. Health index score calculated for each sample.
#' @export freq_health_index
#' @author Rotem Hadar
freq_health_index <- function(exp, ns){
  if(missing(ns)) ns <- import_nsf()
  # filter only relevant ASV, and pivot_longer
  exp <- exp %>% data.frame() %>% 
    tibble::rownames_to_column("Feature_ID") %>% 
    dplyr::right_join(ns) %>% 
    tidyr::pivot_longer(contains("DB"), values_to = "reads", names_to = "sampleid") 
  # calculate HI
  res <- exp %>% 
    dplyr::group_by(sampleid, dir) %>% dplyr::summarise(freq = sum(reads, na.rm = T)) %>% 
    tidyr::pivot_wider(id_cols = sampleid, names_from = dir, values_from = freq) %>% 
    dplyr::transmute(freq_HI = log2((down+.01) / (up+.01)))
  res
}


#' Calculate the ranked health index for each sample in an experiment.
#' @param exp data.frame of ASV exeriment. ASV are rows and samples are the columns.
#' @param ns data.frame of "good" and "bad" bacteria with "dir".
#' @param thresh floa. Threshold to consider level of bacteria absent/not.
#' @return data.frame. Binary table with score for each sample.
#' @export binary_health_index
ranked_health_index <- function(exp, ns, thresh = 0){
  res <- binary_health_index(exp, ns, thresh)
  res <- data.frame(sampleid = res$sampleid, ranked_HI = rank(res$freq_HI))
}

#' Convert AmpliconExperiment type from python calour to standard data.frame ASV table.
#' @param AmpliconExperiment calour object to transorm to data.frame.
#' @author Rotem Hadar
#' @return data.frame features are rows. samples are columns.
#' @export calour_to_df
calour_to_df <- function(AmpliconExperiment){
  df <- data.frame(AmpliconExperiment$data)
  rownames(df) <-  allbiom_exp$sample_metadata$`_sample_id`
  colnames(df) <-  AmpliconExperiment$feature_metadata$`_feature_id`
  data.frame(t(df))
}


