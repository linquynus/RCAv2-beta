#' Compute pairwise DE genes for supervised clustering result.
#'
#' @param rca.obj RCA object.
#' @param logFoldChange Log fold change required to call gene DE.
#' @param method Denotes which test to use. Available options are:
#' \itemize{
#'  \item{"wilcox"} : Identifies differentially expressed genes between two
#'  groups of cells using a Wilcoxon Rank Sum test (default)
#'  \item{"bimod"} : Likelihood-ratio test for single cell gene expression,
#'  (McDavid et al., Bioinformatics, 2013)
#'  \item{"roc"} : Identifies 'markers' of gene expression using ROC analysis.
#'  For each gene, evaluates (using AUC) a classifier built on that gene alone,
#'  to classify between two groups of cells. An AUC value of 1 means that
#'  expression values for this gene alone can perfectly classify the two
#'  groupings (i.e. Each of the cells in cells.1 exhibit a higher level than
#'  each of the cells in cells.2). An AUC value of 0 also means there is perfect
#'  classification, but in the other direction. A value of 0.5 implies that
#'  the gene has no predictive power to classify the two groups. Returns a
#'  'predictive power' (abs(AUC-0.5) * 2) ranked matrix of putative differentially
#'  expressed genes.
#'  \item{"t"} : Identify differentially expressed genes between two groups of
#'  cells using the Student's t-test.
#'  }
#' @param mean.Exp Minimum mean expression of a gene to be considered in the DE gene calculation
#' @param deepsplit If hclust was used for clustering, the desired deepsplit can be specified here.. Values can range from 0 to 4. Default is 1.
#' @param min.pct  only test genes that are detected in a minimum fraction of
#' min.pct cells in either of the two populations. Meant to speed up the function
#' by not testing genes that are very infrequently expressed. Default is 0.25
#' @param min.diff.pct  only test genes that show a minimum difference in the
#' fraction of detection between the two groups. Set to -Inf by default
#' @param random.seed Random seed for downsampling. default is 1
#' @param min.cells.group Minimum number of cells in one of the groups
#' @param pseudocount.use Pseudocount to add to averaged expression values when calculating logFC. 1 by default.
#' @param p.adjust.methods correction method for calculating qvalue. default is BH (or FDR)
#' @param top.genes.per.cluster Number of top DE genes to be considered per cluster
#' @param pairwise Flag indicating whether DE genes should be compared derived in pairwise manner or 1 cluster vs all others (Default).

#' @return RCA object.
#' @export
#'
dataDE <- function(rca.obj,
                   logFoldChange = 1.5,
                   method = "wilcox",
                   mean.Exp = 0.5,
                   deepsplit = 1,
                   min.pct = 0.25,
                   min.diff.pct = -Inf,
                   random.seed = 1,
                   min.cells.group = 3,
                   pseudocount.use = 1,
                   p.adjust.methods =  "BH",
                   top.genes.per.cluster = 10,
                   pairwise=FALSE) {
    df <- c()
    temp.exp = expm1(x = rca.obj$data)
    temp.exp.row = Matrix::rowMeans(temp.exp)
    temp.exp.row = sort(temp.exp.row, decreasing = T)
    temp.exp.row = temp.exp.row[6:length(temp.exp.row)]
    MeanExprsThrs = mean(temp.exp.row)
    ############################
    #hclust used for clustering#
    ############################
    if (class(rca.obj$clustering.out$cellTree) == "hclust") {
        clusters <- rca.obj$clustering.out$dynamicColorsList[[deepsplit]]
        total.clus <- length(unique(clusters))
        remap <- c(1:total.clus)
        names(remap) <-
            unique(rca.obj$clustering.out$dynamicColorsList[[deepsplit]])
        clusters <- as.numeric(as.character(remap[clusters]))
    } else{
        #############################
        #Graph based clustering used#
        #############################
        clusters <- rca.obj$clustering.out$dynamicColorsList[[1]]
        total.clus <- length(unique(clusters))
        remap <- c(1:total.clus)
        names(remap) <- unique(rca.obj$clustering.out$dynamicColorsList[[1]])
        clusters <- as.numeric(as.character(remap[clusters]))
    }
    ###########################
    #Compute pairwise DE genes#
    ###########################
    if (pairwise){
        for (clusteri in 1:(total.clus - 1)) {
            for (clusterj in (clusteri + 1):(total.clus)) {
                print(c(clusteri, clusterj))
                cells.1 <- colnames(rca.obj$data)[which(clusters == clusteri)]
                cells.2 <-
                    colnames(rca.obj$data)[which(clusters == clusterj)]
                marker.genes = ComputePairWiseDE(
                    object = rca.obj$data,
                    cells.1 = cells.1,
                    cells.2 = cells.2,
                    features = NULL,
                    logfc.threshold = logFoldChange,
                    test.use = method,
                    min.pct = min.pct,
                    min.diff.pct = min.diff.pct,
                    verbose = TRUE,
                    only.pos = FALSE,
                    max.cells.per.ident = Inf,
                    random.seed = random.seed,
                    min.cells.group = min.cells.group,
                    pseudocount.use = pseudocount.use,
                    MeanExprsThrs = MeanExprsThrs,
                    p.adjust.methods = p.adjust.methods)
                if (!(is.null(marker.genes))) {
                    if (colnames(marker.genes)[1] != 'myAUC') {
                        marker.genes = marker.genes[marker.genes$p_val_adj < 0.05, ]
                    }
                    if(nrow(marker.genes) > 0) {
                        marker.genes$group1 = clusteri
                        marker.genes$group2 = clusterj
                        marker.genes$gene = rownames(marker.genes)
                        df = rbind(df, marker.genes)
                    }
                }
            }
        }
    }else{
        for (clusteri in 1:(total.clus)) {
            print(clusteri)
            cells.1 <- colnames(rca.obj$data)[which(clusters == clusteri)]
            cells.2 <- colnames(rca.obj$data)[which(clusters != clusteri)]
            marker.genes = ComputePairWiseDE(
                object = rca.obj$data,
                cells.1 = cells.1,
                cells.2 = cells.2,
                features = NULL,
                logfc.threshold = logFoldChange,
                test.use = method,
                min.pct = min.pct,
                min.diff.pct = min.diff.pct,
                verbose = TRUE,
                only.pos = FALSE,
                max.cells.per.ident = Inf,
                random.seed = random.seed,
                min.cells.group = min.cells.group,
                pseudocount.use = pseudocount.use,
                MeanExprsThrs = MeanExprsThrs,
                p.adjust.methods = p.adjust.methods
            )
            if (!(is.null(marker.genes))) {
                if (colnames(marker.genes)[1] != 'myAUC') {
                    marker.genes = marker.genes[marker.genes$p_val_adj < 0.05, ]
                }
                if(nrow(marker.genes) > 0) {
                    marker.genes$group1 = clusteri
                    marker.genes$gene = rownames(marker.genes)
                    df = rbind(df, marker.genes)
                }
            }
        }
    }
    ######################################
    #Determine top x DE genes per Cluster#
    ######################################
    if(pairwise){
        mC1 <- df %>% group_by(group1, group2) %>% dplyr::top_n(n = top.genes.per.cluster, wt = avg_logFC)
        markers2 <- df
        markers2$avg_logFC <- (-1) * (markers2$avg_logFC)
        mC2 <- markers2 %>% group_by(group2, group1) %>% dplyr::top_n(n = top.genes.per.cluster, wt = avg_logFC)
        topMarkers <- data.frame(rbind(
            cbind(Cluster = mC1$group1, Gene = mC1$gene),
            cbind(Cluster = mC2$group2, Gene = mC2$gene)
        ))
        topMarkers <- topMarkers %>% group_by(Cluster) %>% distinct(.keep_all = T)
        topMarkers <- topMarkers[order(topMarkers$Cluster, decreasing = F), ]
    }else{
        topMarkers <- data.frame(Cluster = df$group1, Gene = df$gene, avg_logFC = df$avg_logFC)
        topMarkers <- topMarkers %>% group_by(Cluster) %>% dplyr::top_n(n = top.genes.per.cluster, wt = avg_logFC) %>% distinct(.keep_all = T)
        topMarkers <- topMarkers[order(topMarkers$Cluster, decreasing = F), ]
    }
    rca.obj$DE.genes <- list(All.DE.genes = df, Top.DE.genes = topMarkers)
    return(rca.obj)
}

