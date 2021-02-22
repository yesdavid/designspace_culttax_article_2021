library(cluster)
library(ggplot2)
library(readr)
library(Momocs)
library(dplyr)
library(maptree)
library(ggimage)
library(phangorn)
library(rworldmap)
library(raster)
library(NbClust)

if (!requireNamespace("ggtree", quietly = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
    BiocManager::install("ggtree")
    library(ggtree)
  } else {
    BiocManager::install("ggtree")
    library(ggtree)
  }
} else {library(ggtree)}


rm(list=ls())


outlines_combined_nicolas_2016 <- readRDS(file ="./1_data/outlines_combined_nicholas_2016.RDS")


outlines_combined_nicolas_2016_centered <- Momocs::coo_centre(outlines_combined_nicolas_2016) # center
outlines_combined_nicolas_2016_centered_scaled <- Momocs::coo_scale(outlines_combined_nicolas_2016_centered) # scale
outlines_combined_nicolas_2016_centered_scaled <- Momocs::coo_slidedirection(outlines_combined_nicolas_2016_centered_scaled, 
                                                                             direction = "up")
# Momocs::panel(outlines_combined_nicolas_2016_centered_scaled)


# unification of outlines with catalogue-dataframe
nicolas_fleches_2016_catalog_ids_coordinates <- readr::read_csv(file = "./1_data/nicolas_fleches_2016_catalog_ids_with_coordinates.csv")

outlines_combined_nicolas_2016_names <- names(outlines_combined_nicolas_2016_centered_scaled)
outlines_combined_nicolas_2016_names_splitted <- strsplit(outlines_combined_nicolas_2016_names, split = "_")

ID_and_artefact_ID_list <- list()
for (name_index in 1:length(outlines_combined_nicolas_2016_names)){
  
  ID_and_artefact_ID_interrim_df <- data.frame(ID = paste0(outlines_combined_nicolas_2016_names_splitted[[name_index]][1], "-", outlines_combined_nicolas_2016_names_splitted[[name_index]][2]),
                                               ID_artefact <- outlines_combined_nicolas_2016_names[[name_index]])
  names(ID_and_artefact_ID_interrim_df) <- c("ID", "ID_artefact")
  ID_and_artefact_ID_list[[name_index]] <- ID_and_artefact_ID_interrim_df
  
}

ID_and_artefact_ID_df <- do.call("rbind", ID_and_artefact_ID_list)

nicolas_fleches_2016_catalog_ids_coordinates_artefact_ID <- dplyr::inner_join(ID_and_artefact_ID_df, nicolas_fleches_2016_catalog_ids_coordinates, by = "ID")


outlines_combined_nicolas_2016_centered_scaled <- Momocs::Out(outlines_combined_nicolas_2016_centered_scaled$coo, 
                                                              fac = nicolas_fleches_2016_catalog_ids_coordinates_artefact_ID)
outlines_combined_nicolas_2016_centered_scaled <- Momocs::filter(outlines_combined_nicolas_2016_centered_scaled, 
                                                                 !ID_artefact %in% c("UK_60_XX_pseudo_no_10", "UK_15_XX_pseudo_no_4")) # remove fragmented outliers

length(which(nicolas_fleches_2016_catalog_ids_coordinates_artefact_ID$country == "Denmark")) # number of artefacts from denmark
length(which(nicolas_fleches_2016_catalog_ids_coordinates_artefact_ID$country == "France")) # number of artefacts from france
length(which(nicolas_fleches_2016_catalog_ids_coordinates_artefact_ID$country == "United Kingdom")) # number of artefacts from the uk

# harmonic calibration
outlines_combined_nicolas_2016_centered_scaled_harmonics <- Momocs::calibrate_harmonicpower_efourier(outlines_combined_nicolas_2016_centered_scaled, 
                                                                                                     plot = F)  # Estimates the number of harmonics required for the Fourier methods implemented in Momocs.
outlines_combined_nicolas_2016_centered_scaled_harmonics

# efourier
outlines_combined_nicolas_2016_centered_scaled_efourier <- Momocs::efourier(outlines_combined_nicolas_2016_centered_scaled,
                                                                            nb.h = as.matrix(outlines_combined_nicolas_2016_centered_scaled_harmonics[["minh"]])[[4,1]], # harmonics for 99.9%
                                                                            norm = F) # you selected `norm=TRUE`, which is not recommended. See ?efourier --> probably no problem in our case


# PCA
outlines_combined_nicolas_2016_centered_scaled_PCA <- Momocs::PCA(outlines_combined_nicolas_2016_centered_scaled_efourier) # PCA on Coe objects, using prcomp.

minimum_no_of_pcs_nicolas <- Momocs::scree_min(outlines_combined_nicolas_2016_centered_scaled_PCA,
                                               prop = 0.95) # minimum number of axis to use to retain a given proportion (i.e. prop = 0.99 to describe 99% of the variation) -- reduces computing time in the phylogeny estimation step:
Momocs::scree_plot(outlines_combined_nicolas_2016_centered_scaled_PCA)

# nice PCA plot 
# Create groups
pch.group <- c(rep(21, times=length(which(outlines_combined_nicolas_2016_centered_scaled_PCA$fac$country == "Denmark"))), 
               rep(22, times=length(which(outlines_combined_nicolas_2016_centered_scaled_PCA$fac$country == "France"))),
               rep(23, times=length(which(outlines_combined_nicolas_2016_centered_scaled_PCA$fac$country == "United Kingdom"))))
col.group <- c(rep("skyblue2", times=length(which(outlines_combined_nicolas_2016_centered_scaled_PCA$fac$country == "Denmark"))), 
               rep("gold", times=length(which(outlines_combined_nicolas_2016_centered_scaled_PCA$fac$country == "France"))),
               rep("green2", times=length(which(outlines_combined_nicolas_2016_centered_scaled_PCA$fac$country == "United Kingdom"))))

# plot
svg("./3_output/late_neolithic_early_bronze_age_arrowheads/nicolas_2016_pca_w_outlier_names.svg",
    width = 8,
    height = 8)
plot(outlines_combined_nicolas_2016_centered_scaled_PCA$x[,1],
     outlines_combined_nicolas_2016_centered_scaled_PCA$x[,2],
     xlab=paste("PCA 1 (", round(summary(outlines_combined_nicolas_2016_centered_scaled_PCA)$importance[2]*100, 1), "%)", sep = ""),
     ylab=paste("PCA 2 (", round(summary(outlines_combined_nicolas_2016_centered_scaled_PCA)$importance[5]*100, 1), "%)", sep = ""),
     pch=pch.group,
     col="black",
     bg=col.group,
     cex=1,
     las=1,
     asp=1,
     main = "Nicolas 2016; PCA of Bell Beaker arrowhead outline shapes")

# Add grid lines
abline(v=0, lty=2, col="grey50")
abline(h=0, lty=2, col="grey50")

# Add legend
legend("bottomleft", 
       legend=c("Denmark", "France", "United Kingdom"), 
       col="black", 
       pt.bg=c("skyblue2", "gold", "green2"), 
       pch=c(21, 22, 24), 
       pt.cex=1.5)
dev.off()


#### remove outliers
nicolas_outliers_db <- fpc::dbscan(outlines_combined_nicolas_2016_centered_scaled_PCA$x, eps = 0.3, MinPts = 3)
plot(nicolas_outliers_db, outlines_combined_nicolas_2016_centered_scaled_PCA$x, main = "DBSCAN", frame = FALSE)

nicolas_outliers_cluster <- data.frame(name = row.names(outlines_combined_nicolas_2016_centered_scaled_PCA$x), 
                                       value = nicolas_outliers_db$cluster, 
                                       row.names = NULL)

nicolas_outliers_cluster_outlier_names <- subset(nicolas_outliers_cluster, value != 1)

outlines_combined_nicolas_2016_centered_scaled_PCA$fac$outlier_names <- NA
outlines_combined_nicolas_2016_centered_scaled$fac$cluster <- as.factor(nicolas_outliers_cluster$value)

for (vector_index in 1:length(match(nicolas_outliers_cluster_outlier_names$name, outlines_combined_nicolas_2016_centered_scaled_PCA$fac$ID_artefact))){
  
  current_index <- match(nicolas_outliers_cluster_outlier_names$name, outlines_combined_nicolas_2016_centered_scaled_PCA$fac$ID_artefact)[vector_index]
  
  outlines_combined_nicolas_2016_centered_scaled_PCA$fac$outlier_names[current_index] <- nicolas_outliers_cluster_outlier_names$name[vector_index]
}


# look at the shape of the outliers
nicolas_2016_with_outliers <- Momocs::slice(outlines_combined_nicolas_2016_centered_scaled, match(nicolas_outliers_cluster_outlier_names$name, outlines_combined_nicolas_2016_centered_scaled_PCA$fac$ID_artefact))

Momocs::panel(nicolas_2016_with_outliers,
              fac = "country_code",
              names = T,
              col = "grey")

# analysis without outliers

nicolas_2016_without_outliers <- Momocs::slice(outlines_combined_nicolas_2016_centered_scaled, -match(nicolas_outliers_cluster_outlier_names$name, outlines_combined_nicolas_2016_centered_scaled_PCA$fac$ID_artefact))


# harmonic calibration
nicolas_2016_without_outliers_harmonics <- Momocs::calibrate_harmonicpower_efourier(nicolas_2016_without_outliers, 
                                                                                    plot = F)  # Estimates the number of harmonics required for the Fourier methods implemented in Momocs. This is the only step in this section that produces data we need in the subsequent step.

# efourier
nicolas_2016_without_outliers_efourier <- Momocs::efourier(nicolas_2016_without_outliers,
                                                           nb.h = as.matrix(nicolas_2016_without_outliers_harmonics[["minh"]])[[4,1]], # harmonics for 99.9%
                                                           norm = F) 
# PCA
nicolas_2016_without_outliers_PCA <- Momocs::PCA(nicolas_2016_without_outliers_efourier) # PCA on Coe objects, using prcomp.

nicolas_2016_screeplot_wo_outliers <- Momocs::scree_plot(nicolas_2016_without_outliers_PCA, nax = 1:8)
ggsave(filename = "./3_output/late_neolithic_early_bronze_age_arrowheads/nicolas_2016_screeplot_wo_outliers.svg",
       nicolas_2016_screeplot_wo_outliers,
       width = 20,
       height = 20,
       units = "cm")
ggsave(filename = "./3_output/late_neolithic_early_bronze_age_arrowheads/nicolas_2016_screeplot_wo_outliers.png",
       nicolas_2016_screeplot_wo_outliers,
       width = 20,
       height = 20,
       units = "cm")

minimum_no_of_pcs_nicolas_without_outliers <- ncol(nicolas_2016_without_outliers_PCA$x) #use all PC axes


gg_nicolas <- Momocs::PCcontrib(nicolas_2016_without_outliers_PCA,
                                nax = 1:5,
                                sd.r = c(-2,-1,0,1,2))
library(ggplot2)
gg_nicolas <- gg_nicolas$gg +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

ggsave(filename = "./3_output/late_neolithic_early_bronze_age_arrowheads/nicolas_2016_pca_variance.svg",
       gg_nicolas,
       width = 20,
       height = 20,
       units = "cm")
ggsave(filename = "./3_output/late_neolithic_early_bronze_age_arrowheads/nicolas_2016_pca_variance.png",
       gg_nicolas,
       width = 20,
       height = 20,
       units = "cm")

# nice PCA plot 
# Create groups
pch.group <- c(rep(21, times=length(which(nicolas_2016_without_outliers_PCA$fac$country == "Denmark"))), 
               rep(22, times=length(which(nicolas_2016_without_outliers_PCA$fac$country == "France"))),
               rep(23, times=length(which(nicolas_2016_without_outliers_PCA$fac$country == "United Kingdom"))))
col.group <- c(rep("skyblue2", times=length(which(nicolas_2016_without_outliers_PCA$fac$country == "Denmark"))), 
               rep("gold", times=length(which(nicolas_2016_without_outliers_PCA$fac$country == "France"))),
               rep("green2", times=length(which(nicolas_2016_without_outliers_PCA$fac$country == "United Kingdom"))))

svg("./3_output/late_neolithic_early_bronze_age_arrowheads/nicolas_2016_pca_without_outlier_names.svg",
    width = 8, height = 8)
# plot
plot(nicolas_2016_without_outliers_PCA$x[,1],
     nicolas_2016_without_outliers_PCA$x[,2],
     xlab=paste("PCA 1 (", round(summary(nicolas_2016_without_outliers_PCA)$importance[2]*100, 1), "%)", sep = ""),
     ylab=paste("PCA 2 (", round(summary(nicolas_2016_without_outliers_PCA)$importance[5]*100, 1), "%)", sep = ""),
     pch=pch.group,
     col="black",
     bg=col.group,
     cex=1,
     las=1,
     asp=1,
     main = "Nicolas 2016; PCA of Bell Beaker arrowhead outline shapes without outliers")
# Add grid lines
abline(v=0, lty=2, col="grey50")
abline(h=0, lty=2, col="grey50")
# Add legend
legend("bottomleft", 
       legend=c("Denmark", "France", "United Kingdom"), 
       col="black", 
       pt.bg=c("skyblue2", "gold", "green2"), 
       pch=c(21, 22, 24), 
       pt.cex=1.5)
dev.off()



##################################### Clustering ##################################### 


## distance matrix
nicolas_2016_without_outliers_PCA_as_matrix <- as.matrix(nicolas_2016_without_outliers_PCA$x[,1:minimum_no_of_pcs_nicolas_without_outliers])
nicolas_2016_without_outliers_dist <- dist(nicolas_2016_without_outliers_PCA_as_matrix)

## hierarchical clustering
nicolas_2016_without_outliers_upgma <- hclust(nicolas_2016_without_outliers_dist,
                                              method = "average")
nicolas_2016_without_outliers_upgmc <- hclust(nicolas_2016_without_outliers_dist,
                                              method = "centroid")
nicolas_2016_without_outliers_wardD2 <- hclust(nicolas_2016_without_outliers_dist,
                                               method = "ward.D2")

nicolas_list_different_clustering_methods <- list("UPGMA" = nicolas_2016_without_outliers_upgma,
                                                  "UPGMC" = nicolas_2016_without_outliers_upgmc,
                                                  "Ward" = nicolas_2016_without_outliers_wardD2)


############# k 
library(NbClust)
nicolas_2016_without_outliers_NbClust_ward <- NbClust::NbClust(data = nicolas_2016_without_outliers_PCA$x[,1:minimum_no_of_pcs_nicolas_without_outliers],
                                                               distance = "euclidean",
                                                               method = "ward.D2",
                                                               index = c("gap", "silhouette"))
nicolas_2016_without_outliers_NbClust <- nicolas_2016_without_outliers_NbClust_ward$All.index
nicolas_2016_without_outliers_NbClust_df <- as.data.frame(nicolas_2016_without_outliers_NbClust)
nicolas_2016_without_outliers_NbClust_df$NClust <- 1:nrow(nicolas_2016_without_outliers_NbClust)+1
silhouette_plot <- ggplot(nicolas_2016_without_outliers_NbClust_df, aes(x = NClust, y = Silhouette)) + geom_point() + geom_line() + theme_bw() +
  scale_x_continuous(breaks = seq(2, max(nicolas_2016_without_outliers_NbClust_df$NClust)+1, by = 1)) +
  xlab("Number of clusters") +
  ylab("Average silhouette value")
gap_plot <- ggplot(nicolas_2016_without_outliers_NbClust_df, aes(x = NClust, y = Gap)) + geom_point() + geom_line() + theme_bw() +
  scale_x_continuous(breaks = seq(2, max(nicolas_2016_without_outliers_NbClust_df$NClust)+1, by = 1)) +
  xlab("Number of clusters") +
  ylab("Gap statistic (k)")
cowplot::plot_grid(silhouette_plot, gap_plot)

height_silhouette_plot <- 100
width_silhouette_plot <- 100

ggsave(silhouette_plot,
       filename = "./3_output/late_neolithic_early_bronze_age_arrowheads/Ward/silhouette_plot_NbClust_wardD2.svg",
       height  = height_silhouette_plot,
       width = width_silhouette_plot,
       units = "mm")
ggsave(silhouette_plot,
       filename = "./3_output/late_neolithic_early_bronze_age_arrowheads/Ward/silhouette_plot_NbClust_wardD2.png",
       height  = height_silhouette_plot,
       width = width_silhouette_plot,
       units = "mm")

n_clusters_nicolas_2016_without_outliers <- 9

set.seed(1)
nicolas_colors_without_outliers <- randomcoloR::distinctColorPalette(n_clusters_nicolas_2016_without_outliers)


#### create pruned dendrogram with mean shapes, artefact panel per cluster and spatial distribution map of created clusters for list of fusion methods established above (UPGMA, UPGMC, Ward's method)


world <- rgeos::gBuffer(rworldmap::getMap(resolution = "high"), byid=TRUE, width=0)
clipper_europe <- as(raster::extent(-7, 13, 46, 58), "SpatialPolygons")
proj4string(clipper_europe) <- CRS(proj4string(world))
world_clip <- raster::intersect(world, clipper_europe)
world_clip_f <- fortify(world_clip)


layout = "circular"
img_size <- 0.02
tip_lab_size <- 3

nicolas_treecuts_by_algorithm_list <- list()

for (current_algorithm_name in names(nicolas_list_different_clustering_methods)){
  
  current_algorithms_clustering <- nicolas_list_different_clustering_methods[[current_algorithm_name]]
  
  
  
  ############### tree without outliers
  
  current_treecut <- data.frame(ID_artefact = row.names(as.data.frame(cutree(current_algorithms_clustering,
                                                                             k = n_clusters_nicolas_2016_without_outliers))),
                                cluster = as.factor(as.data.frame(cutree(current_algorithms_clustering, 
                                                                         k = n_clusters_nicolas_2016_without_outliers))[[1]]))
  nicolas_treecuts_by_algorithm_list[[current_algorithm_name]] <- current_treecut
  
  nicolas_2016_without_outliers_w_cluster <- Momocs::Out(nicolas_2016_without_outliers$coo,
                                                         fac = current_treecut)
  # harmonic calibration
  nicolas_2016_without_outliers_w_cluster_harmonics <- Momocs::calibrate_harmonicpower_efourier(nicolas_2016_without_outliers_w_cluster, 
                                                                                                plot = F)
  # efourier
  nicolas_2016_without_outliers_w_cluster_efourier <- Momocs::efourier(nicolas_2016_without_outliers_w_cluster,
                                                                       nb.h = as.matrix(nicolas_2016_without_outliers_w_cluster_harmonics[["minh"]])[[4,1]], # harmonics for 99.9%
                                                                       norm = F) 
  # PCA
  nicolas_2016_without_outliers_w_cluster_PCA <- Momocs::PCA(nicolas_2016_without_outliers_w_cluster_efourier) # PCA on Coe objects, using prcomp.
  
  
  
  
  ##############
  # panel plots by cluster without outliers
  ##############
  
  nicolas_algorithm_folder <- paste0("./3_output/late_neolithic_early_bronze_age_arrowheads/", current_algorithm_name, "/")
  dir.create(nicolas_algorithm_folder)
  
  ape::write.nexus(ape::as.phylo(current_algorithms_clustering),
                   file = paste0(nicolas_algorithm_folder, current_algorithm_name, "_full_tree.tre"))
  
  nicolas_algorithm_pathname <- paste0(nicolas_algorithm_folder, current_algorithm_name, "_clusters_k", n_clusters_nicolas_2016_without_outliers,"_no_outliers/")
  nicolas_algorithm_pathname_panel <- paste0(nicolas_algorithm_pathname, "panels/")
  nicolas_algorithm_pathname_means <- paste0(nicolas_algorithm_pathname, "means/")
  
  dir.create(path = nicolas_algorithm_pathname)
  dir.create(path = nicolas_algorithm_pathname_panel)
  dir.create(path = nicolas_algorithm_pathname_means)
  
  # PCA colored by current algorithm cluster
  png(file=paste0(nicolas_algorithm_pathname, "pairs_PCA_w_k", n_clusters_nicolas_2016_without_outliers),
      width = 800, height = 600, units = "px")
  pairs(nicolas_2016_without_outliers_w_cluster_PCA$x[,c(1:5)],
        col = nicolas_colors_without_outliers[nicolas_2016_without_outliers_w_cluster_PCA$fac$cluster],
        cex = 1,
        pch = 15,
        lower.panel = NULL)
  dev.off()
  
  svg(file=paste0(nicolas_algorithm_pathname, "PCA_w_k", n_clusters_nicolas_2016_without_outliers, ".svg"),
      width = 8, height = 6)
  plot(nicolas_2016_without_outliers_w_cluster_PCA$x[,1],
       nicolas_2016_without_outliers_w_cluster_PCA$x[,2],
       xlab=paste("PCA 1 (", round(summary(nicolas_2016_without_outliers_w_cluster_PCA)$importance[2]*100, 0), "%)", sep = ""),
       ylab=paste("PCA 2 (", round(summary(nicolas_2016_without_outliers_w_cluster_PCA)$importance[5]*100, 0), "%)", sep = ""),
       bg=nicolas_colors_without_outliers[nicolas_2016_without_outliers_w_cluster_PCA$fac$cluster], 
       pch = 21,
       cex=1,
       las=1,
       asp=1)
  # Add grid lines
  abline(v=0, lty=2, col="grey50")
  abline(h=0, lty=2, col="grey50")
  dev.off()
  
  # panels for each cluster
  for (i in 1:n_clusters_nicolas_2016_without_outliers){
    
    mypath_png <- paste0(nicolas_algorithm_pathname_panel, "nicolas_2016_without_outliers_w_cluster_colors_cluster_", i, ".png")
    
    png(file=mypath_png,
        width = 800, height = 800, units = "px")
    
    Momocs::panel(Momocs::slice(nicolas_2016_without_outliers_w_cluster, cluster == i),
                  main = paste("Cluster", i),
                  col = nicolas_colors_without_outliers[i])
    
    dev.off()
    
    mypath_svg <- paste0(nicolas_algorithm_pathname_panel, "nicolas_2016_without_outliers_w_cluster_colors_cluster_", i, ".svg")
    
    svg(file=mypath_svg, width = 8, height = 8)
    
    Momocs::panel(Momocs::slice(nicolas_2016_without_outliers_w_cluster, cluster == i),
                  main = paste("Cluster", i),
                  col = nicolas_colors_without_outliers[i])
    
    dev.off()
  }
  
  
  
  # mean shapes
  min_no_of_coordinates <- list()
  for (cluster_index in 1:n_clusters_nicolas_2016_without_outliers){
    
    current_shapes <- Momocs::slice(nicolas_2016_without_outliers_w_cluster, cluster == cluster_index)
    
    min_no_of_coordinates[[cluster_index]] <- rep(NA, length(current_shapes))
    
    for (i in 1:length(current_shapes)){
      min_no_of_coordinates[[cluster_index]][i] <- nrow(current_shapes$coo[[i]])
    }
    
    min_no_of_coordinates[[cluster_index]] <- min(min_no_of_coordinates[[cluster_index]])
    
  }
  
  ## shapes get INTERPOLATED to common number of landmarks (lowest number of landmarks per cluster)
  mean_shapes_cluster_without_outliers <- list()
  for (cluster_index in 1:n_clusters_nicolas_2016_without_outliers){
    mean_shapes_cluster_without_outliers[[cluster_index]] <- Momocs::MSHAPES(Momocs::coo_interpolate(Momocs::slice(nicolas_2016_without_outliers_w_cluster, cluster == cluster_index),
                                                                                                     n = min_no_of_coordinates[[cluster_index]])$coo)
  }
  
  nicolas_2016_without_outliers_w_cluster_PCA_mean_shapes_cluster_out <- Momocs::Out(mean_shapes_cluster_without_outliers,
                                                                                     fac = data.frame(cluster = paste0("cluster_", c(1:n_clusters_nicolas_2016_without_outliers))))
  
  
  
  for (i in 1:n_clusters_nicolas_2016_without_outliers){
    
    mypath_png <- paste0(nicolas_algorithm_pathname_means, paste0("Cluster ", i, " (n=",colSums(table(current_treecut))[[i]], ")"), "_mean_shp.png")
    
    png(file=mypath_png,
        width = 800, height = 800, units = "px")
    
    Momocs::panel(Momocs::slice(nicolas_2016_without_outliers_w_cluster_PCA_mean_shapes_cluster_out, cluster == paste0("cluster_",i)),
                  main = NULL,
                  col = nicolas_colors_without_outliers[i])
    
    dev.off()
    
    
    mypath_svg <- paste0(nicolas_algorithm_pathname_means, paste0("Cluster ", i, " (n=",colSums(table(current_treecut))[[i]], ")"), "_mean_shp.svg")
    
    svg(file=mypath_svg, width = 8, height = 8)
    
    Momocs::panel(Momocs::slice(nicolas_2016_without_outliers_w_cluster_PCA_mean_shapes_cluster_out, cluster == paste0("cluster_",i)),
                  main = NULL,
                  col = nicolas_colors_without_outliers[i])
    
    dev.off()
  }
  
  if (current_algorithm_name != "UPGMC"){
    
    # pruned tree with mean shapes
    current_algorithms_clustering_cuttree <- maptree::clip.clust(current_algorithms_clustering, 
                                                                 current_treecut, 
                                                                 k=n_clusters_nicolas_2016_without_outliers)
    ape::write.nexus(ape::as.phylo(current_algorithms_clustering_cuttree),
                     file = paste0(nicolas_algorithm_folder, current_algorithm_name, "_pruned_tree.tre"))
    
    ######################## change tip labels of chosen trees
    new_labels <- rep(NA, length(current_algorithms_clustering_cuttree$labels))
    for (i in 1:length(current_algorithms_clustering_cuttree$labels)){
      new_labels[i] <- paste0("Cluster ", current_algorithms_clustering_cuttree$labels[i], " (n=",colSums(table(current_treecut))[[i]], ")")
    }
    
    
    current_algorithms_clustering_cuttree$labels <- new_labels
    
    
    # 
    current_pruned_tree <- ggtree(current_algorithms_clustering_cuttree, 
                                  layout = "rectangular") + 
      xlim(NA,6) +
      geom_tiplab(aes(image=paste0(nicolas_algorithm_pathname_means, label, "_mean_shp.png")), 
                  geom="image", offset=0.5, align=2, hjust = 0.5) + 
      geom_tiplab(geom='label', offset=1.5, hjust=.5, size = 5) + 
      geom_treescale()
    
    ggsave(filename = paste0(nicolas_algorithm_pathname, "pruned_dendro.svg"),
           plot = current_pruned_tree,
           device = "svg",
           width = 25,
           height = 35,
           units = "cm")
  }
  
  
  
  ### distribution map
  data_wo_outliers <- nicolas_2016_without_outliers$fac
  data_wo_outliers$cluster <- NULL
  
  data_wo_outliers_unique_sites <- dplyr::distinct(data_wo_outliers, ID, .keep_all =T)
  
  data_wo_outliers <- dplyr::left_join(data_wo_outliers, nicolas_2016_without_outliers_w_cluster$fac, by = "ID_artefact")
  
  
  no_of_clusters <- length(unique(data_wo_outliers$cluster))
  cluster_names <- c()
  
  for (i in 1:no_of_clusters){
    cluster_names[i] <- paste0("Cluster ", i)
    
  }
  names(cluster_names) <- 1:no_of_clusters
  
  
  current_distribution_map <- ggplot() +
    geom_polygon(data = world_clip_f, 
                 aes(x = long, y = lat, group = group),
                 fill = NA, colour = "grey") + 
    geom_point(data = data_wo_outliers_unique_sites[,c("site", "lng", "lat")],  
               aes(x = lng, y = lat, alpha = 0.7), 
               shape = 3, color = "black") +
    geom_jitter(data = data_wo_outliers,
                aes(x = lng, y = lat,
                    fill = cluster), 
                shape = 21, size = 3,
                width = 0.1, height = 0.1
    ) +   
    facet_wrap(~cluster,
               scales = "fixed", 
               labeller = as_labeller(cluster_names)) +
    coord_quickmap() +  
    theme_classic() + 
    xlab("Longitude") +
    ylab("Latitude") +
    theme(legend.position = "none")
  
  ggsave(filename = paste0(nicolas_algorithm_pathname, "distribution_map_", current_algorithm_name, ".svg"),
         plot = current_distribution_map,
         device = "svg",
         width = 30,
         height = 25,
         units = "cm")
  ggsave(filename = paste0(nicolas_algorithm_pathname, "distribution_map_", current_algorithm_name, ".png"),
         plot = current_distribution_map,
         device = "png",
         width = 30,
         height = 25,
         units = "cm")
}





############################# NJ typochronology ############################# 


typochronologie_csv <- readr::read_csv("./1_data/nicolas_2017_typochronologie.csv")

typochronologie_csv <- dplyr::distinct(typochronologie_csv, ID_country, .keep_all = T)

typochronologie_csv_UK <- subset(typochronologie_csv, Country == "UK")
typochronologie_csv_FR <- subset(typochronologie_csv, Country == "FR")


nicolas_2016_without_outliers_PCA_as_df <- as.data.frame(nicolas_2016_without_outliers_PCA$x[,1:minimum_no_of_pcs_nicolas_without_outliers])


nicolas_2016_without_outliers_PCA_as_df$ID_country <- sapply(rownames(nicolas_2016_without_outliers_PCA_as_df),
                                                             function(x){
                                                               strsplit(x, split = "(?<=.{5})", perl = TRUE)[[1]][1] # RegEx "(?<=.{5})" means: split into chunks of 5 characters long
                                                             })  

########## FR
nicolas_2016_without_outliers_PCA_as_df_subset_typochron_FR <- subset(nicolas_2016_without_outliers_PCA_as_df, ID_country %in% typochronologie_csv_FR$ID_country)

names_artefacts_ID <- data.frame(artefact_ID = rownames(nicolas_2016_without_outliers_PCA_as_df_subset_typochron_FR),
                                 ID_country = nicolas_2016_without_outliers_PCA_as_df_subset_typochron_FR$ID_country)
names_artefacts_ID_and_period <- dplyr::left_join(names_artefacts_ID, typochronologie_csv_FR[,c("ID_country", "Period")], by = "ID_country")
names_artefacts_ID_and_period$Period <- as.factor(names_artefacts_ID_and_period$Period)

nicolas_2016_without_outliers_PCA_as_df_subset_typochron_FR$ID_country <- NULL
nicolas_typochron_FR_dist <- dist(nicolas_2016_without_outliers_PCA_as_df_subset_typochron_FR, method = "euclidean")

nicolas_typochron_FR_NJ <- phangorn::NJ(nicolas_typochron_FR_dist)

additional_information_period <- data.frame(Period = names_artefacts_ID_and_period[,c("Period")])
rownames(additional_information_period) <- names_artefacts_ID_and_period$artefact_ID



NJ_nicolas_ggtree <- ggtree(nicolas_typochron_FR_NJ) %<+% names_artefacts_ID_and_period +
  geom_tiplab(size=3, 
              aes(label = ID_country,
                  color = ID_country)) +
  geom_treescale() + 
  scale_colour_discrete(na.translate = F) + 
  guides(color=guide_legend(title="Site ID"))


NJ_nicolas_gheatmap <- gheatmap(NJ_nicolas_ggtree,
                                additional_information_period,
                                offset=0.1,
                                width=0.1,
                                legend_title="Period",
                                colnames = T) +
  scale_x_ggtree() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  guides(fill=guide_legend(title="Period")) +
  scale_fill_discrete(na.translate = F)

NJ_nicolas_gheatmap_height <- 10
NJ_nicolas_gheatmap_width <- 10

ggsave(NJ_nicolas_gheatmap,
       filename = "./3_output/late_neolithic_early_bronze_age_arrowheads/NJ_with_heatmap_periods.svg",
       width = NJ_nicolas_gheatmap_width,
       height = NJ_nicolas_gheatmap_height,
       device = "svg",
       units = "in")
ggsave(NJ_nicolas_gheatmap,
       filename = "./3_output/late_neolithic_early_bronze_age_arrowheads/NJ_with_heatmap_periods.png",
       width = NJ_nicolas_gheatmap_width,
       height = NJ_nicolas_gheatmap_height,
       device = "png",
       units = "in")
