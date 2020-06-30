#' Perform cluster specific QC
#'
#' @param rca.obj RCA object
#' @param cluster.labels vector of cluster labels
#' @param width width of plot in inches. Default is 20.
#' @param height height of plot in inches. Default is 20.
#' @param folderpath path to save heatmap to
#' @param filename file name of saved heatmap
#'
#' @export
#'

plotClusterQuality <- function(rca.obj, cluster.labels, width, height, folderpath = ".", filename = "RCA_Cluster_Quality.pdf") {

    if(missing(width)){
        width <- 9
    }

    if(missing(height)){
        height <- 3*length(unique(cluster.labels))
    }

    # Create data frame for cell-cluster mapping
    cluster.df <- data.frame(Cell = colnames(rca.obj$data), Cluster = cluster.labels)

    # Create empty data frame for QC parameters
    quality.df <- data.frame(Cluster = character(), nGene = numeric(), nUMI = numeric(), pMito = numeric(), stringsAsFactors = FALSE)

    # For each cluster
    for(cluster in unique(cluster.labels)) {

        # Subset cluster data
        data <- rca.obj$raw.data[, subset(cluster.df$Cell, cluster.df$Cluster == cluster), drop = FALSE]

        # Compute nGene vector
        nGeneVec <- Matrix::colSums(data>0)

        # Compute nUMI vector
        nUMIVec <- Matrix::colSums(data)

        # Select mito genes
        mito.genes = grep(pattern = "^MT-", x = rownames(data), value = T)

        # Compute percent.mito vector
        pMitoVec <- Matrix::colSums(data[mito.genes, , drop = FALSE])/Matrix::colSums(data)

        # Append QC vectors to data frame
        cluster.quality.df <- data.frame(Cluster = rep(cluster, ncol(data)), nGene = nGeneVec, nUMI = nUMIVec, pMito = pMitoVec)
        quality.df <- rbind(quality.df, cluster.quality.df)
    }



    ## added by Quy
    # plot list for nGene vs nUMI
    plot_list1 <- list()
    # plot list for nGene vs pMito
    plot_list2 <- list()
    # plot list for nUMI vs pMito
    plot_list3 <- list()
    for (i in 1:length(unique(cluster.labels))){
        cluster_i <- unique(cluster.labels)[i]
        quality.df.i <- quality.df[which(quality.df$Cluster == cluster_i),, drop=FALSE]
        if (nrow(quality.df.i) >1) {
            p.i1 <- ggplot2::ggplot(data = quality.df.i, ggplot2::aes(x = nGene, y = nUMI)) +
                ggplot2::geom_point(size = 1) + ggplot2::geom_jitter() +
                ggplot2::theme_bw() + ggplot2::ggtitle(cluster_i)+
                ggplot2::stat_density_2d()


            p.i2 <- ggplot2::ggplot(data = quality.df.i, ggplot2::aes(x = nGene, y = pMito)) +
                ggplot2::geom_point(size = 1) + ggplot2::geom_jitter() +
                ggplot2::theme_bw() + ggplot2::ggtitle(cluster_i)+
                ggplot2::stat_density_2d()


            p.i3 <- ggplot2::ggplot(data = quality.df.i, ggplot2::aes(x = nUMI, y = pMito)) +
                ggplot2::geom_point(size = 1) + ggplot2::geom_jitter() +
                ggplot2::theme_bw() + ggplot2::ggtitle(cluster_i)+
                ggplot2::stat_density_2d()

        } else {
            p.i1 <- ggplot2::ggplot(data = quality.df.i, ggplot2::aes(x = nGene, y = nUMI)) +
                ggplot2::geom_point(size = 1) +
                ggplot2::theme_bw() + ggplot2::ggtitle(cluster_i)


            p.i2 <- ggplot2::ggplot(data = quality.df.i, ggplot2::aes(x = nGene, y = pMito)) +
                ggplot2::geom_point(size = 1) +
                ggplot2::theme_bw() + ggplot2::ggtitle(cluster_i)


            p.i3 <- ggplot2::ggplot(data = quality.df.i, ggplot2::aes(x = nUMI, y = pMito)) +
                ggplot2::geom_point(size = 1) +
                ggplot2::theme_bw() + ggplot2::ggtitle(cluster_i)
        }

        plot_list1[[i]] <- p.i1
        plot_list2[[i]] <- p.i2
        plot_list3[[i]] <- p.i3
    }


    p1 <- arrangeGrob(grobs=plot_list1,
                      nrow=length(unique(cluster.labels)))
    p2 <- arrangeGrob(grobs=plot_list2,
                      nrow=length(unique(cluster.labels)))
    p3 <- arrangeGrob(grobs=plot_list3,
                      nrow=length(unique(cluster.labels)))

    pdf(file = paste0(folderpath, "/", "Cluster_quality_", filename), width = width, height = height)
    grid.arrange(p1, p2, p3, ncol=3)
    dev.off()

}
