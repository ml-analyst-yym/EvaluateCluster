#' @title Evaluate cluster effectiveness
#' @description Help choose the best clustering resolution in Seurat.
#' @param data A SeuratObject containing multiple resolution clusters
#' @param cores Windows default 1. other system is free
#' @param rogue_span The parameter of rogue controls the loess. Default 0.6.Recommend 0.6-0.8
#' @returns Return a list include the result of all methods
#' @import Seurat
#' @import dplyr
#' @importFrom clusterCrit intCriteria
#' @import parallel
#' @import ROGUE
#' @import lisi
#' @import mrtree
#' @examples
#' \dontrun{
#' EvaluateClusters(data, cores=40)
#' EvaluateClusters(data, cores=8, rogue_span=0.8)
#' EvaluateClusters(data, cores=1, rogue_span=0.5)
#' }
#' @export
EvaluateCluster <- function(data, cores = 1, rogue_span=0.6) {
  set.seed(12315)
  
  # ---- Data preparation ----
  data_ori <- as.data.frame(data@reductions$umap@cell.embeddings)
  meta_ori <- data@meta.data
  label <- grep("^RNA_snn_res\\.", colnames(meta_ori), value = TRUE)
  
  
  # ---- ROGUE ----

  message("ROGUE start \n")
  expr <- GetAssayData(data, assay = "RNA", layer = "counts")
  rogue_list <- parallel::mclapply(label, function(res_col) {
    ROGUE::rogue(expr, platform = "UMI",
                 labels = meta_ori[[res_col]],
                 samples = meta_ori$orig.ident,
                 span = rogue_span)
  }, mc.cores = cores)
  
  names(rogue_list) <- label
  rogue_res <- sapply(rogue_list, function(df) mean(as.matrix(df), na.rm = TRUE))
  names(rogue_res) <- sub("RNA_snn_res\\.", "", names(rogue_res))
  message("ROGUE finish \n")
  
  
  # ---- LISI ----
 
  message("LISI start \n")
  lisi_i <- compute_lisi(data_ori, meta_ori, label)
  lisi_res <- apply(lisi_i,2,mean)
  lisi_res <- 1/lisi_res
  names(lisi_res) <- sub("RNA_snn_res\\.", "", label)
  message("LISI finish \n")
  
  
  # ---- MRtree ----
  
  message("MRtree start \n")
  out=mrtree(data,prefix = "RNA_snn_res.",n.cores=cores)
  stab.out = stability_plot(out)
  ks.flat = apply(out$labelmat.flat, 2, FUN=function(x) length(unique(x)))
  matched_cols <- names(ks.flat)[ks.flat %in% stab.out$df$resolution]
  matched_scores <- stab.out$df$ari[match(ks.flat[matched_cols], stab.out$df$resolution)]
  mrtree_r <- data.frame(
    colname = matched_cols,
    resolution = ks.flat[matched_cols],
    ari = matched_scores)
  mrtree_r[,1] <- sub("RNA_snn_res\\.", "", mrtree_r[,1])
  mrtree_res <- mrtree_r[,3]
  names(mrtree_res) <- sub("RNA_snn_res\\.", "", label)
  message("MRtree finish \n")
  
  # ---- Internal indices ----
  
  message("Internal indices start \n")
  metrics <- purrr::map_dfr(label, function(res_col) {
    cluster_labels <- meta_ori[[res_col]]
    valid <- !is.na(cluster_labels)
    cluster_labels <- cluster_labels[valid]
    if (length(unique(cluster_labels)) <= 1) {
      return(data.frame(res = sub("RNA_snn_res\\.", "", res_col),
                        CH = NA, Sil = NA, DB = NA))
    }
    crit <- clusterCrit::intCriteria(as.matrix(data_ori[valid,]),
                                     as.integer(as.factor(cluster_labels)),
                                     c("Calinski_Harabasz", "Silhouette", "Davies_Bouldin"))
    data.frame(res = sub("RNA_snn_res\\.", "", res_col),
               CH = crit$calinski_harabasz,
               Sil = crit$silhouette,
               DB = crit$davies_bouldin)
  })
  
  CH_res <- metrics$CH; names(CH_res) <- metrics$res
  Sil_res <- metrics$Sil; names(Sil_res) <- metrics$res
  DB_res <- metrics$DB; names(DB_res) <- metrics$res
  message("Internal indices finish \n")
  
  # ---- Rank aggregation ----
  rank_CH  <- rank(-CH_res, na.last = "keep")
  rank_Sil <- rank(-Sil_res, na.last = "keep")
  rank_DB  <- rank(DB_res, na.last = "keep")
  rank_avg <- rowMeans(cbind(rank_CH, rank_Sil, rank_DB), na.rm = TRUE)
  names(rank_avg) <- metrics$res
  
  
  # ---- Combine results ----
  
  result_df <- list(
    Resolution = as.numeric(metrics$res),
    ROGUE = rogue_res,
    LISI = lisi_res,
    MRtree = mrtree_res,
    CH = CH_res,
    Silhouette = Sil_res,
    DB = DB_res,
    Rank = rank_avg
  )
  
  return(list(scores = result_df,
              meta_ori = meta_ori))
}




#' @title Evaluate cluster effectiveness by AUC
#' @description Calculate the AUC of Markers to help choose the max resolution.
#' @param data A SeuratObject containing multiple resolution clusters.
#' @param cores Windows default 1. other system is free.
#' @param auc_cutoff The threshold of AUC score, default 0.6. It will affect the maximum resolution.
#' @param min_pct The parameter of FindAllMarkers, default 0.05. It will affect the result of AUC.Not recommended to change
#' @param logfc_threshold The parameter of FindAllMarkers, default 0.1. It will affect the result of AUC.Not recommended to change
#' @returns Return a list include the result of all methods
#' @import Seurat
#' @import dplyr
#' @import parallel
#' @examples
#' \dontrun{
#' CalculateAUC(data, cores=40)
#' CalculateAUC(data, cores=8, min_pct = 0.05, logfc_threshold = 0.2)
#' CalculateAUC(data, cores=1, min_pct = 0.1, logfc_threshold = 0.1)
#' }
#' @export
CalculateAUC <- function(data, cores = 1, min_pct = 0.05, logfc_threshold = 0.1, auc_cutoff = 0.6){
  
  message("AUC start \n")
  res_cols <- grep("^RNA_snn_res\\.", colnames(data@meta.data), value = TRUE)
  
  cl <- parallel::makeCluster(cores)
  parallel::clusterExport(cl, varlist = c("data","res_cols","min_pct","logfc_threshold","auc_cutoff"),envir = environment())
  parallel::clusterEvalQ(cl, {
    library(Seurat)
    library(dplyr)
  })
  
  marker_counts_list <- parallel::parLapply(cl, res_cols, function(res_col){
    scRNA_tmp <- data
    scRNA_tmp <- SetIdent(scRNA_tmp, value = res_col)
    markers <- FindAllMarkers(
      scRNA_tmp,
      test.use = "roc",
      min.pct = min_pct,
      logfc.threshold = logfc_threshold
    ) %>%
      filter(myAUC > auc_cutoff) %>%
      count(cluster, name = "number") %>%
      mutate(resolution = res_col)
    rm(scRNA_tmp)
    gc()
    return(markers)
  })
  parallel::stopCluster(cl)
  
  marker_counts <- dplyr::bind_rows(marker_counts_list)
  message("AUC finish \n")
  return(marker_counts)
}




#' @title Choose the resolution
#' @description Choose the best resolution based on the result of other function.
#' @param result_eva the result of EvaluateCluster
#' @param marker_counts the result of CalculateAUC
#' @param rogue_threshold The threshold of rogue score, default 0.8. It will affect the minimum resolution
#' @returns Return a list include the result of all methods
#' @import purrr
#' @examples
#' \dontrun{
#' Choose_res(result_eva, marker_counts)
#' Choose_res(result_eva, marker_counts, rogue_threshold = 0.8)
#' Choose_res(result_eva, marker_counts, rogue_threshold = 0.9)
#' }
#' @export
Choose_res <- function(result_eva, marker_counts, rogue_threshold = 0.85){
  
  result_df = result_eva$scores
  rogue_res = result_df$ROGUE
  lisi_res = result_df$LISI
  mrtree_res = result_df$MRtree
  CH_res = result_df$CH
  Sil_res = result_df$Silhouette
  DB_res = result_df$DB
  rank_avg=  result_df$Rank
  res_cols <- grep("^RNA_snn_res\\.", colnames(result_eva$meta_ori), value = TRUE)
  
  # The lower limit resolution
  rogue_valid <- names(rogue_res[rogue_res > rogue_threshold])
  min_rogue_res <- if (length(rogue_valid) > 0) min(as.numeric(rogue_valid)) else NA
  
  # The upper limit resolution (markerless cluster detected)
  res_numeric <- sort(result_df$Resolution)
  auc_upper_res <- {
    idx <- purrr::map_lgl(res_cols, function(res) {
      clusters <- sort(as.numeric(as.character(unique(result_eva$meta_ori[[res]]))))
      detected <- sort(as.numeric(as.character(unique(marker_counts$cluster[marker_counts$resolution == res]))))
      any(!clusters %in% detected)
    })
    if (any(idx)) {
      current_res <- as.numeric(sub("RNA_snn_res\\.", "", res_cols[which(idx)[1]]))
      res_index <- match(current_res, res_numeric)
      if (!is.na(res_index) && res_index > 1) {
        res_numeric[res_index - 1]
      } else if (res_index == 1) {
        current_res
      }else NA
    } else {max(res_numeric)}
  }
  
  # LISI + MRtree mean, select top3
  mean_lisi_mrtree <- (lisi_res + mrtree_res)/2
  names(mean_lisi_mrtree) <- names(lisi_res)
  
  all_res_names <- names(mean_lisi_mrtree)
  allowed_mask <- rep(TRUE, length(all_res_names))
  # apply ROGUE lower bound if present
  if (!is.na(min_rogue_res)) {
    allowed_mask <- allowed_mask & (as.numeric(all_res_names) >= as.numeric(min_rogue_res))}
  # apply AUC upper bound if present
  if (!is.na(auc_upper_res)) {
    allowed_mask <- allowed_mask & (as.numeric(all_res_names) <= as.numeric(auc_upper_res))}
  allowed_res <- all_res_names[allowed_mask]
  
  # if filtering yields none, fall back to unfiltered set
  if (length(allowed_res) == 0) {
    consider_vec <- mean_lisi_mrtree
  } else {consider_vec <- mean_lisi_mrtree[allowed_res]}
  
  # pick top5 from considered vector
  ord <- order(consider_vec, decreasing = TRUE, na.last = NA)
  if (length(ord) == 0) {
    top5_candidates <- character(0)
  } else {
    topN <- min(5, length(ord))
    top5_candidates <- names(consider_vec)[ord[1:topN]]
  }
  
  # final selection among top3 by internal rank_avg
  final_rank <- rank_avg[names(rank_avg) %in% top5_candidates]
  if (length(final_rank) == 0) {
    # fallback: if no topN candidates (unlikely), pick best by rank across allowed_res if available, else overall
    if (length(allowed_res) > 0) {
      best_res <- allowed_res[which.min(rank_avg[allowed_res])]
    } else {best_res <- names(which.min(rank_avg))}
  } else {best_res <- names(which.min(final_rank))}
  

  
  return(list(scores = result_df,
              auc_markers = marker_counts,
              rogue_min_res = min_rogue_res,
              auc_upper_res = auc_upper_res,
              best_resolution = best_res))
}


  

