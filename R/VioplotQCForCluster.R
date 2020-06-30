#' Fviolin plot of QC metrics
#'
#' @param rca.obj RCA object.
#' @param prefix_name prefix filename.
#' @export
#'
VioplotQCForCluster <- function(rca.obj, prefix_name="1"){
    unique.label <- unique(rca.obj$clustering.out$dynamicColorsList[[1]])
    for (i in 1:length(unique.label)){
        label.i <- unique.label[i]
        cluster.df <- data.frame(Cell = colnames(rca.obj$data),
                                 Cluster = rca.obj$clustering.out$dynamicColorsList[[1]])
        subset.cluster <- rca.obj$raw.data[, subset(cluster.df$Cell, cluster.df$Cluster == label.i),
                                           drop = FALSE]
        # Compute nGene vector
        nGeneVec.i <- Matrix::colSums(subset.cluster>0)

        # Compute nUMI vector
        nUMIVec.i <- Matrix::colSums(subset.cluster)

        # Select mito genes
        mito.genes.i = grep(pattern = "^MT-", x = rownames(subset.cluster),
                            value = T)

        # Compute percent.mito vector
        pMitoVec.i <- Matrix::colSums(subset.cluster[mito.genes.i, , drop = FALSE])/Matrix::colSums(subset.cluster)

        cluster.quality.df.i <- data.frame(Cluster = rep(label.i, ncol(subset.cluster)),
                                           nGene = nGeneVec.i, nUMI = nUMIVec.i,
                                           pMito = pMitoVec.i)
        cluster.quality.df.i$pos <- 1
        g1 <- ggplot(data = cluster.quality.df.i, aes(x=pos,y = nGeneVec.i))+geom_violin(trim = FALSE)+
            geom_jitter() + xlab("nGene")+ylab("counts")+ggtitle(label.i)
        g2 <- ggplot(data = cluster.quality.df.i, aes(x=pos,y = nUMIVec.i))+geom_violin(trim = FALSE)+
            geom_jitter() + xlab("nUMI")+ylab("counts")+ggtitle(label.i)
        g3 <- ggplot(data = cluster.quality.df.i, aes(x=pos,y = pMitoVec.i))+geom_violin(trim = FALSE)+
            geom_jitter() + xlab("pMt")+ylab("counts")+ggtitle(label.i)

        pdf(file = paste0(prefix_name,"_cluster_quality_for_",i,
                          "_cluster_as_color_",label.i,".pdf"),
            width = 6, height = 4)
        grid.arrange(g1, g2, g3, nrow = 1)
        dev.off()
    }
}
