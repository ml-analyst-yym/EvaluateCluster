#' @title Evaluate cluster effectiveness
#' @description Help choose the best clustering resolution in Seurat.
#' @param data A SeuratObject containing multiple resolution clusters
#' @param cores Windows default 1. other system is free
#' @param rogue First method, default false, recommend Windows to False
#' @param lisi Second method, default false
#' @param mrtree Third method, default false, recommend Windows to False
#' @param CH_index Fourth method, default true
#' @param plot Determine whether the function gives the plot of result
#'
#' @returns Return a list include the result of 4 methods
#' @export
#' @importFrom SeuratObject JoinLayers
#' @importFrom SeuratObject GetAssayData
#' @import tidyverse
#' @importFrom graphics abline
#' @importFrom parallel mclapply
#' @examples
#' \dontrun{
#' EvaluateClusters(data,cores=20)
#' EvaluateClusters(data,cores=8,rogue=FALSE,plot=FALSE)
#' EvaluateClusters(data,rogue=FALSE,mrtree=FALSE,lisi=FALSE,cores=1)
#' }
EvaluateCluster <- function(data, cores = 1,
                             rogue = FALSE, lisi = FALSE, mrtree = FALSE, CH_index = TRUE,
                             plot = TRUE) {
  set.seed(12315)

  # ---- Data preparation ----
  data = JoinLayers(data)
  data_ori <- as.data.frame(data@reductions$umap@cell.embeddings)
  meta_ori <- data@meta.data
  label <- grep("^RNA_snn_res\\.", colnames(meta_ori), value = TRUE)

  # ---- ROGUE ----
  if (rogue) {
    if (!requireNamespace("ROGUE", quietly = TRUE)) {
      stop("ROGUE package is required to run this method. Please install from GitHub.")
    }
  }
  if (rogue) {
    cat("ROGUE start \n")
    expr <- GetAssayData(data, assay = "RNA", layer = "counts")
    if (cores > 1) {
      rogue_list <- parallel::mclapply(label, function(res_col) {
        ROGUE::rogue(expr, platform = "UMI",
              labels = meta_ori[[res_col]],
              samples = meta_ori$orig.ident,
              span = 0.6)
      }, mc.cores = cores)
    } else {
      rogue_list <- lapply(label, function(res_col) {
        ROGUE::rogue(expr, platform = "UMI",
              labels = meta_ori[[res_col]],
              samples = meta_ori$orig.ident,
              span = 0.6)
      })
    }
    names(rogue_list) <- label
    rogue_res <- sapply(rogue_list, function(df) mean(as.matrix(df), na.rm = TRUE))
    names(rogue_res) <- sub("RNA_snn_res\\.", "", names(rogue_res))
    cat("ROGUE finish \n")
  } else {rogue_res <- NULL}

  # ---- LISI ----
  if (lisi) {
    if (!requireNamespace("lisi", quietly = TRUE)) {
      stop("lisi package is required to run this method. Please install from GitHub.")
    }
    # 调用 ROGUE::rogue(...)
  }
  if (lisi) {
    cat("LISI start \n")
    lisi_i <- lisi::compute_lisi(data_ori, meta_ori, label)
    lisi_res <- apply(lisi_i, 2, mean)
    names(lisi_res) <- sub("RNA_snn_res\\.", "", names(lisi_res))
    lisi_res <- 1/lisi_res
    cat("LISI finish \n")
  } else {lisi_res <- NULL}

  # ---- MRtree ----
  if (mrtree) {
    if (!requireNamespace("mrtree", quietly = TRUE)) {
      stop("mrtree package is required to run this method. Please install from GitHub.")
    }
    # 调用 ROGUE::rogue(...)
  }
  if (mrtree) {
    cat("MRtree start \n")
    out = mrtree::mrtree(data, prefix = "RNA_snn_res.", n.cores = cores)
    stab.out = mrtree::stability_plot(out)
    ks.flat = apply(out$labelmat.flat, 2, FUN = function(x) length(unique(x)))
    matched_cols <- names(ks.flat)[ks.flat %in% stab.out$df$resolution]
    matched_scores <- stab.out$df$ari[match(ks.flat[matched_cols], stab.out$df$resolution)]
    mrtree_r <- data.frame(
      colname = matched_cols,
      resolution = ks.flat[matched_cols],
      ari = matched_scores)
    mrtree_r[,1] <- sub("RNA_snn_res\\.", "", mrtree_r[,1])
    mrtree_res <- mrtree_r[,3]
    names(mrtree_res) <- mrtree_r[,1]
    cat("MRtree finish \n")
  } else {mrtree_res <- NULL}

  # ---- CH_index ----
  if (CH_index) {
    cat("CH_index start \n")
    CH_i <- function(df_coord, cluster) {
      df <- cbind(df_coord, cluster = cluster)
      n <- nrow(df)
      k <- length(unique(cluster))
      if (k == 1 | k == n) return(NA)
      overall_centroid <- colMeans(df_coord)
      Wk <- sum(sapply(unique(cluster), function(cl) {
        members <- df[df$cluster == cl, 1:ncol(df_coord)]
        cl_centroid <- colMeans(members)
        sum(rowSums((members - matrix(cl_centroid, nrow(members), ncol(df_coord), byrow = TRUE))^2))
      }))
      Bk <- sum(sapply(unique(cluster), function(cl) {
        members <- df[df$cluster == cl, 1:ncol(df_coord)]
        cl_centroid <- colMeans(members)
        n_cl <- nrow(members)
        n_cl * sum((cl_centroid - overall_centroid)^2)
      }))
      CH <- (Bk/(k-1))/(Wk/(n-k))
      CH_avg <- CH/n
      return(CH_avg)
    }
    df_cluster <- meta_ori[, label, drop = FALSE]
    CH_res <- sapply(df_cluster, function(cl) CH_i(data_ori, cl))
    names(CH_res) <- sub("RNA_snn_res\\.", "", names(CH_res))
    cat("CH_index finish \n")
  } else {CH_res <- NULL}

  # ---- result ----
  result_list <- list()
  if (!is.null(CH_res))    result_list$CH_index <- CH_res
  if (!is.null(lisi_res))  result_list$LISI     <- lisi_res
  if (!is.null(rogue_res)) result_list$ROGUE    <- rogue_res
  if (!is.null(mrtree_res))result_list$MRtree   <- mrtree_res

  if (length(result_list) == 0) {
    warning("No score calculated.")
    return(NULL)
  }

  result_df <- do.call(rbind, lapply(names(result_list), function(nm) {
    x <- result_list[[nm]]
    all_names <- unique(unlist(lapply(result_list, names)))
    vals <- rep(NA_real_, length(all_names))
    names(vals) <- all_names
    vals[names(x)] <- as.numeric(x)
    rowdf <- as.data.frame(t(vals), stringsAsFactors = FALSE)
    rownames(rowdf) <- nm
    rowdf
  }))
  result_df <- as.data.frame(result_df)
  comb <- list(combined_score = colMeans(result_df,na.rm = T))
  cat("Score:\n")
  print(comb$combined_score)
  cat("\nBest reduction: ", names(which.max(comb$combined_score)), "\n")

  # ---- plot ----
  if (plot) {
    plot(as.numeric(names(comb$combined_score)),
         comb$combined_score, type = "b", pch = 19,
         xlab = "Resolution", ylab = "Combined Score",
         main = "Overall score")
    abline(v = as.numeric(names(which.max(comb$combined_score))), col = "red", lty = 2)
  }

  return(list(score = result_df,
              combined_score = comb$combined_score,
              best_resolution = names(which.max(comb$combined_score))))
}
































