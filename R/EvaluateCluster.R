#' @title Evaluate cluster effectiveness
#' @description Help choose the best clustering resolution in Seurat.
#' @param data A SeuratObject containing multiple resolution clusters
#' @param cores Windows default 1. other system is free
#' @param rogue_threshold The threshold of rogue score, default 0.8. It will affect the minimum resolution
#' @param auc_cutoff The threshold of AUC score, default 0.6. It will affect the maximum resolution
#'
#' @returns Return a list include the result of all methods
#' @export
#' @import Seurat
#' @import purrr
#' @import dplyr
#' @importFrom clusterCrit intCriteria
#' @import parallel
#' @import ROGUE
#' @import lisi
#' @import mrtree
#' @examples
#' \dontrun{
#' EvaluateClusters(data, cores=20)
#' EvaluateClusters(data, cores=8, rogue_threshold=0.8, auc_cutoff=0.6)
#' EvaluateClusters(data, cores=1, rogue_threshold=0.85, auc_cutoff=0.7)
#' }

EvaluateCluster <- function(data, cores = 1, auc_cutoff = 0.6,
                            rogue_threshold = 0.8) {
  set.seed(12315)
  
  # ---- Data preparation ----
  data <- JoinLayers(data)
  data_ori <- as.data.frame(data@reductions$umap@cell.embeddings)
  meta_ori <- data@meta.data
  label <- grep("^RNA_snn_res\\.", colnames(meta_ori), value = TRUE)
  
  
  # ---- ROGUE ----
  if (!requireNamespace("ROGUE", quietly = TRUE)) {
    message("Package 'ROGUE' not found. You can install it with:")
    message("  devtools::install_github('PaulingLiu/ROGUE')")
    stop("Please install the ROGUE package before running EvaluateCluster().")
  }
  message("ROGUE start \n")
  expr <- GetAssayData(data, assay = "RNA", layer = "counts")
  rogue_list <- parallel::mclapply(label, function(res_col) {
    ROGUE::rogue(expr, platform = "UMI",
                 labels = meta_ori[[res_col]],
                 samples = meta_ori$orig.ident,
                 span = 0.6)
  }, mc.cores = cores)
  
  names(rogue_list) <- label
  rogue_res <- sapply(rogue_list, function(df) mean(as.matrix(df), na.rm = TRUE))
  names(rogue_res) <- sub("RNA_snn_res\\.", "", names(rogue_res))
  message("ROGUE finish \n")
  
  rogue_valid <- names(rogue_res[rogue_res > rogue_threshold])
  min_rogue_res <- if (length(rogue_valid) > 0) min(as.numeric(rogue_valid)) else NA
  
  
  # ---- AUC ----
  message("AUC start \n")
  res_cols <- label
  
  cl <- parallel::makeCluster(cores)
  parallel::clusterExport(cl, varlist = c("data", "res_cols", "auc_cutoff"),envir = environment())
  parallel::clusterEvalQ(cl, {
    library(Seurat)
    library(dplyr)
  })
  
  marker_counts_list <- parallel::parLapply(cl, res_cols, function(res_col){
    scRNA_tmp <- data
    scRNA_tmp <- SetIdent(scRNA_tmp, value = res_col)
    markers <- FindAllMarkers(
      scRNA_tmp,
      test.use = "roc"
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
  
  # ---- Identify upper limit resolution (markerless cluster detected) ----
  res_numeric <- sort(as.numeric(sub("RNA_snn_res\\.", "", res_cols)))
  auc_upper_res <- {
    idx <- purrr::map_lgl(res_cols, function(res) {
      clusters <- sort(as.numeric(as.character(unique(meta_ori[[res]]))))
      detected <- sort(as.numeric(as.character(unique(marker_counts$cluster[marker_counts$resolution == res]))))
      !identical(clusters, seq(0, max(clusters)))
    })
    if (any(idx)) {
      current_res <- as.numeric(sub("RNA_snn_res\\.", "", res_cols[which(idx)[1]]))
      res_index <- match(current_res, res_numeric)
      if (!is.na(res_index) && res_index > 1) {
        res_numeric[res_index - 1]
      } else if (res_index == 1) {
        current_res
      }else NA
    } else NA
  }
  message("AUC finish \n")
  
  # ---- LISI ----
  if (!requireNamespace("lisi", quietly = TRUE)) {
    message("Package 'lisi' not found. You can install it with:")
    message("  devtools::install_github('immunogenomics/lisi')")
    stop("Please install the lisi package before running EvaluateCluster().")
  }
  message("LISI start \n")
  lisi_res <- sapply(label, function(res_col) {
    lisi_mat <- lisi::compute_lisi(data_ori, meta_ori, res_col)
    mean(lisi_mat[, res_col], na.rm = TRUE)
  })
  names(lisi_res) <- sub("RNA_snn_res\\.", "", label)
  message("LISI finish \n")
  
  # ---- MRtree ----
  if (!requireNamespace("mrtree", quietly = TRUE)) {
    message("Package 'mrtree' not found. You can install it with:")
    message("  devtools::install_github('pengminshi/mrtree')")
    stop("Please install the mrtree package before running EvaluateCluster().")
  }
  message("MRtree start \n")
  mrtree_res <- sapply(label, function(res_col) {
    out <- mrtree::mrtree(data, prefix = res_col, n.cores = cores)
    stab.out <- mrtree::stability_plot(out)
    ks.flat <- apply(out$labelmat.flat, 2, FUN = function(x) length(unique(x)))
    matched_cols <- names(ks.flat)[ks.flat %in% stab.out$df$resolution]
    matched_scores <- stab.out$df$ari[match(ks.flat[matched_cols], stab.out$df$resolution)]
    mean(matched_scores, na.rm = TRUE)
  })
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
  names(rank_avg) <- names(CH_res)
  
  # ---- LISI + MRtree mean, select top3 ----
  mean_lisi_mrtree <- mean(lisi_res + mrtree_res)
  top3_candidates <- names(sort(mean_lisi_mrtree, decreasing = TRUE))[1:3]
  final_rank <- rank_avg[names(rank_avg) %in% top3_candidates]
  best_res <- names(which.min(final_rank))
  
  # ---- Combine results ----
  result_df <- data.frame(
    Resolution = as.numeric(names(CH_res)),
    ROGUE = rogue_res[match(names(CH_res), names(rogue_res))],
    LISI = lisi_res[match(names(CH_res), names(lisi_res))],
    MRtree = mrtree_res[match(names(CH_res), names(mrtree_res))],
    CH = CH_res,
    Silhouette = Sil_res,
    DB = DB_res,
    Rank = rank_avg
  )
  
  message("\nBest resolution : ", best_res, "\n")
  
  return(list(
    scores = result_df,
    rogue_min_res = min_rogue_res,
    auc_upper_res = auc_upper_res,
    auc_markers = marker_counts,
    best_resolution = best_res
  ))
}

