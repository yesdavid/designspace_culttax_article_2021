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

output_folder <- file.path(".", "3_output", "final_palaeolithic_earliest_mesolithic_tanged_points")
dir.create(output_folder,
           recursive = T)


tanged_points_tps <- Momocs::import_tps(file.path(".", "1_data", "TPS_TP_27_09_2019.TPS"), 
                                        curves = TRUE)
tanged_points <- readr::read_csv(file.path(".", "1_data", "tanged.points.csv"),
                                 col_types = cols(Context = col_character(),
                                                  Site = col_character()))

outlinetp <- Momocs::Out(tanged_points_tps$coo, 
                         fac = tanged_points)
outlinetp <- Momocs::filter(outlinetp, 
                            Burinated=="N") ### removal of burinated examples

outlinetpTPTC <- Momocs::filter(outlinetp, 
                                ATU %in% "Tanged Point Technocomplex")  
outlinetpATU  <- Momocs::filter(outlinetp, 
                                !ATU %in% "Unspecified")  
outlinetpATU2 <- Momocs::filter(outlinetp, 
                                !ATU2 %in% "Unspecified") 
outlinetpNAT  <- Momocs::filter(outlinetp, 
                                !NAT %in% c("Tanged Point (Unspecified)", NA))
outlinetpNAT_subset  <- Momocs::filter(outlinetpNAT, 
                                       !NAT %in% c("Swiderian Point", "Chwalibogowice Point", "Hamburgian Point", "Teyjat Point")) # we remove these NAT-categories from the dataset, as their sample size is way too small.


##### count NAT_subset occurences
occurences_TP_NAT_subset_list <- list()
occurences_TP_NAT_subset_countries_list <- list()
for (i in 1:length(unique(outlinetpNAT_subset$fac$NAT))){
  current_NAT_subset <- unique(outlinetpNAT_subset$fac$NAT)[i]
  occurences_TP_NAT_subset_list[[i]] <- data.frame(NAT_subset = current_NAT_subset,
                                                   n = length(which(outlinetpNAT_subset$fac$NAT==current_NAT_subset)))
  occurences_TP_NAT_subset_countries_list[[i]] <- unique(subset(outlinetpNAT_subset$fac, NAT==current_NAT_subset)$Country)
}
occurences_TP_NAT_subset <- Reduce(rbind, occurences_TP_NAT_subset_list)
names(occurences_TP_NAT_subset) <- c("NAT", "n")
occurences_TP_NAT_subset

readr::write_csv(occurences_TP_NAT_subset, 
                 file.path(output_folder, "occurences_TP_NAT_subset.csv"))

# unique chronozones
unique(outlinetpNAT_subset$fac$`Chronozone (Relative)`)
unique(outlinetpNAT_subset$fac[,c("NAT", "Chronozone (Relative)")])

## Country X NAT
TP_country_x_nat <- outlinetpNAT_subset$fac[,c("NAT", "Country")]
TP_country_x_nat$NAT <- unlist(strsplit(TP_country_x_nat$NAT, " Point"))
TP_country_x_nat$Country <- factor(TP_country_x_nat$Country, levels = c("France", "Belgium", "The Netherlands", "Germany", "Denmark", "Sweden", "Poland", "Belarus", "Ukraine", "Russia", "Lithuania"))
table_TP_country_x_nat <- t(table(TP_country_x_nat))
print(table_TP_country_x_nat)

## Chronozone X NAT 
TP_chronozone_x_nat <- outlinetpNAT_subset$fac[,c("NAT", "Chronozone (Relative)")]
TP_chronozone_x_nat$NAT <- unlist(strsplit(TP_chronozone_x_nat$NAT, " Point"))
# TP_country_x_nat$`Chronozone (Relative)` <- factor(TP_country_x_nat$`Chronozone (Relative)`, 
#                                                    levels = c())
table_TP_chronozone_x_nat <- t(table(TP_chronozone_x_nat))
print(table_TP_chronozone_x_nat)

#### GMM outline procedure

  outlinetpNAT_subset <- Momocs::coo_centre(outlinetpNAT_subset) # center
  outlinetpNAT_subset <- Momocs::coo_scale(outlinetpNAT_subset) # scale
  outlinetpNAT_subset <- Momocs::coo_slidedirection(outlinetpNAT_subset, direction = "up")

### HARMONIC CALIBRATION ###
  tanged_points_harmonics <- Momocs::calibrate_harmonicpower_efourier(outlinetpNAT_subset, 
                                                                      plot = F)  # Estimates the number of harmonics required for the Fourier methods implemented in Momocs. This is the only step in this section that produces data we need in the subsequent step.
  
  # Momocs::calibrate_reconstructions_efourier(outlinetpNAT_subset) # Calculates and displays reconstructed shapes using a range of harmonic number. Can be used for a visual comparison with the maximal fit.

###  EFOURIER FILE CREATION ###
  tanged_points_efourier <- Momocs::efourier(outlinetpNAT_subset,
                                                           nb.h = as.matrix(tanged_points_harmonics[["minh"]])[[4,1]], # harmonics for 99.9% -> is 99.9% necessary? couldn't it be just 99%?
                                                           norm = F) # see above
########### PCA
  tanged_points_PCA <- Momocs::PCA(tanged_points_efourier) # PCA on Coe objects


  


  
  scree_plot <- Momocs::scree_plot(tanged_points_PCA,
                                   nax = 1:5) +  
    theme_bw() +
    theme(text = element_text(size=20))
  
  minimum_no_of_pcs_tanged_points_PCA <- ncol(tanged_points_PCA$x)
  
  # ggsave(plot = scree_plot,
  #        filename = file.path(output_folder, "PCA_screeplot.svg"),
  #        device = "svg",
  #        width = 20,
  #        height = 15,
  #        units = "cm")
  
  
  gg <- Momocs::PCcontrib(tanged_points_PCA,
                          nax = 1:5,
                          sd.r = c(-2,-1,0,1,2))
  pc_contrib_plot <- gg$gg + 
    theme_bw() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          text = element_text(size=20))
  # ggsave(plot = pc_contrib_plot,
  #        filename = file.path(output_folder, "PCA_contrib.svg"),
  #        device = "svg",
  #        width = 5,
  #        height = 7,
  #        units = "in")

  
  
  ### clustering
  tanged_points_PCA_dist <- dist(tanged_points_PCA$x[,1:minimum_no_of_pcs_tanged_points_PCA])
  
  tanged_points_PCA_upgma <- hclust(tanged_points_PCA_dist,
                                    method = "average")
  
  tanged_points_PCA_upgmc <- hclust(tanged_points_PCA_dist^2, # squared euclidean distance
                                    method = "centroid")
  
  tanged_points_PCA_wardD2 <- hclust(tanged_points_PCA_dist,
                                     method = "ward.D2") # when using ward.D2 distance has not to be squared
  
  chosen_TP_list_different_clustering_methods <- list("UPGMA" = tanged_points_PCA_upgma,
                                                      "UPGMC" = tanged_points_PCA_upgmc,
                                                      "Ward" = tanged_points_PCA_wardD2)

  
  
  
  
  
  # optimal number of _k_ for data using Ward's method
  tanged_points_PCA_NbClust_ward <- NbClust::NbClust(data = tanged_points_PCA$x[,1:minimum_no_of_pcs_tanged_points_PCA],
                                                     distance = "euclidean",
                                                     method = "ward.D2",
                                                     index = c("gap", "silhouette"))
  TP_NbClust <- tanged_points_PCA_NbClust_ward$All.index
  TP_NbClust_df <- as.data.frame(TP_NbClust)
  TP_NbClust_df$NClust <- 1:nrow(TP_NbClust)+1
  silhouette_plot <- ggplot(TP_NbClust_df, aes(x = NClust, y = Silhouette)) + geom_point() + geom_line() + theme_bw() +
    scale_x_continuous(breaks = seq(2, max(TP_NbClust_df$NClust)+1, by = 1)) +
    xlab("Number of clusters") +
    ylab("Average silhouette value") +
    theme(text=element_text(size=20))
  gap_plot <- ggplot(TP_NbClust_df, aes(x = NClust, y = Gap)) + geom_point() + geom_line() + theme_bw() +
    scale_x_continuous(breaks = seq(2, max(TP_NbClust_df$NClust)+1, by = 1)) +
    xlab("Number of clusters") +
    ylab("Gap statistic (k)") +
    theme(text=element_text(size=20))
  cowplot::plot_grid(silhouette_plot, gap_plot)
  
  height_silhouette_plot <- 200
  width_silhouette_plot <- 200
  
  output_folder_ward <- file.path(output_folder, "Ward")
  dir.create(output_folder_ward,
             recursive = T)

  ggsave(silhouette_plot,
         filename = file.path(output_folder_ward, "silhouette_plot_NbClust_wardD2.svg"),
         height  = height_silhouette_plot,
         width = width_silhouette_plot,
         units = "mm")
  ggsave(silhouette_plot,
         filename = file.path(output_folder_ward, "silhouette_plot_NbClust_wardD2.png"),
         height  = height_silhouette_plot,
         width = width_silhouette_plot,
         units = "mm")
  
  
      n_clusters_tanged_points_PCA <- 10 
  paste0("n_clusters_tanged_points_PCA = ", n_clusters_tanged_points_PCA)
  
  set.seed(1)
  tanged_points_PCA_colors <- RColorBrewer::brewer.pal(n = n_clusters_tanged_points_PCA, 
                                                       "Paired")
  
  
  
  
  world <- rgeos::gBuffer(rworldmap::getMap(resolution = "high"), byid=TRUE, width=0)
  clipper_europe <- as(raster::extent(0, 38, 47.5, 58.5), "SpatialPolygons") # -4, 38, 44.5, 57.5
  proj4string(clipper_europe) <- CRS(proj4string(world))
  world_clip <- raster::intersect(world, clipper_europe)
  world_clip_f <- fortify(world_clip)
  
  
  layout = "circular"
  img_size <- 0.02
  tip_lab_size <- 3
  
  
  for (current_algorithm_name in names(chosen_TP_list_different_clustering_methods)){
    
    current_algorithms_clustering <- chosen_TP_list_different_clustering_methods[[current_algorithm_name]]
    
    
    TP_algorithm_folder <- file.path(output_folder, current_algorithm_name)
    dir.create(TP_algorithm_folder)
    TP_algorithm_pathname <- file.path(TP_algorithm_folder, paste0(current_algorithm_name, "_clusters_k", n_clusters_tanged_points_PCA,"_no_outliers"))
    dir.create(path = TP_algorithm_pathname)
    
    
    ############### tree without outliers
    
    current_treecut <- data.frame(ID_artefact = row.names(as.data.frame(cutree(current_algorithms_clustering,
                                                                               k = n_clusters_tanged_points_PCA))),
                                  cluster = as.factor(as.data.frame(cutree(current_algorithms_clustering,
                                                                           k = n_clusters_tanged_points_PCA))[[1]]))
    
    current_treecut <- dplyr::left_join(current_treecut, tanged_points, by = c("ID_artefact" = "FN"))
    
    
    chosen_TP_w_cluster <- Momocs::Out(outlinetpNAT_subset$coo,
                                         fac = current_treecut)

    # harmonic calibration
    chosen_TP_w_cluster_harmonics <- Momocs::calibrate_harmonicpower_efourier(chosen_TP_w_cluster,
                                                                              plot = F)
    # efourier
    chosen_TP_w_cluster_efourier <- Momocs::efourier(chosen_TP_w_cluster,
                                                     nb.h = as.matrix(chosen_TP_w_cluster_harmonics[["minh"]])[[4,1]], # harmonics for 99.9%
                                                     norm = F) 
    # PCA
    chosen_TP_w_cluster_PCA <- Momocs::PCA(chosen_TP_w_cluster_efourier) # PCA on Coe objects, using prcomp.
    

    chosen_TP_w_cluster_PCA_df <- as.data.frame(chosen_TP_w_cluster_PCA$x)
    chosen_TP_w_cluster_PCA_df$cluster <- chosen_TP_w_cluster_PCA$fac$cluster
    
    # chosen_TP_w_cluster_PCA_df <- dplyr::left_join(chosen_TP_w_cluster_PCA_df, chosen_TP_w_cluster_PCA$fac, by = "cluster")
    
    a <- ggplot(data = chosen_TP_w_cluster_PCA_df, aes(fill = cluster)) +
      geom_point(pch = 21, size = 3) +
      coord_fixed(ratio =1) +
      theme_classic() +
      theme(legend.position = "right",
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12),
            legend.text = element_text(size = 14),
            legend.title = element_text(size = 16)) +
      scale_fill_manual(values = tanged_points_PCA_colors) +
      geom_hline(yintercept=0, linetype="dashed", alpha = 0.5) + 
      geom_vline(xintercept=0, linetype="dashed", alpha = 0.5) +
      coord_cartesian(xlim = c(min(chosen_TP_w_cluster_PCA_df[,c("PC1")]) + min(chosen_TP_w_cluster_PCA_df[,c("PC1")])*0.1, 
                               max(chosen_TP_w_cluster_PCA_df[,c("PC1")]) + max(chosen_TP_w_cluster_PCA_df[,c("PC1")])*0.1),
                      ylim = c(min(chosen_TP_w_cluster_PCA_df[,c("PC2")]) + min(chosen_TP_w_cluster_PCA_df[,c("PC2")])*0.1, 
                               max(chosen_TP_w_cluster_PCA_df[,c("PC2")]) + max(chosen_TP_w_cluster_PCA_df[,c("PC2")])*0.1)) +
      theme(text = element_text(size=20))
    
    
    b <- a + aes(x = PC1, y = PC2) +
      xlab(paste0("PC1 (", round(chosen_TP_w_cluster_PCA$eig[1]*100, digits = 0), "%)")) +
      ylab(paste0("PC2 (", round(chosen_TP_w_cluster_PCA$eig[2]*100, digits = 0), "%)")) +
      guides(color = FALSE, fill = FALSE)
    
    c <- a + aes(x = PC1, y = PC3) +
      xlab(paste0("PC1 (", round(chosen_TP_w_cluster_PCA$eig[1]*100, digits = 0), "%)")) +
      ylab(paste0("PC3 (", round(chosen_TP_w_cluster_PCA$eig[3]*100, digits = 0), "%)")) + 
      labs(fill = "cluster")  +
      guides(color = FALSE, fill = FALSE) 
    
    d <- a + aes(x = PC1, y = PC4) +
      xlab(paste0("PC1 (", round(chosen_TP_w_cluster_PCA$eig[1]*100, digits = 0), "%)")) +
      ylab(paste0("PC4 (", round(chosen_TP_w_cluster_PCA$eig[4]*100, digits = 0), "%)")) +
      theme(legend.position="bottom", legend.box = "horizontal") + 
      guides(fill=guide_legend(nrow=2,byrow=TRUE, title = "Cluster"))
    #guides(color = FALSE, fill = FALSE)
    

    
    col_1 <- cowplot::plot_grid(b, c, d,
                       align = "hv", ncol = 1, labels = c("A", "C", "E"), axis = "l")
    col_2 <- cowplot::plot_grid(scree_plot, pc_contrib_plot, rel_heights = c(1,2),
                                align = "hv", ncol = 1, labels = c("B","D"), axis = "l")
    whole_grid <- cowplot::plot_grid(col_1, col_2)
    
    height <- 380
    width <- 285
    ggsave(whole_grid,
           filename = file.path(TP_algorithm_pathname, "whole_grid.svg"),
           height  = height,
           width = width,
           units = "mm")
    ggsave(whole_grid,
           filename = file.path(TP_algorithm_pathname, "whole_grid.png"),
           height  = height,
           width = width,
           units = "mm")
    
    ##############
    # panel plots by cluster without outliers
    ##############
    
    
    TP_algorithm_pathname_panel <- file.path(TP_algorithm_pathname, "panels")
    TP_algorithm_pathname_means <- file.path(TP_algorithm_pathname, "means")
    
    dir.create(path = TP_algorithm_pathname_panel)
    dir.create(path = TP_algorithm_pathname_means)
    
    # PCA pairs plot colored by current algorithm cluster
    png(file=file.path(TP_algorithm_pathname, paste0("pairs_PCA_w_k", n_clusters_tanged_points_PCA, ".png")),
        width = 800, height = 600, units = "px")
    pairs(chosen_TP_w_cluster_PCA$x[,c(1:5)],
          col = tanged_points_PCA_colors[chosen_TP_w_cluster_PCA$fac$cluster],
          cex=1.5, cex.axis = 1.5, cex.lab = 1.5,
          pch = 15,
          lower.panel = NULL)
    dev.off()
    
    svg(file=file.path(TP_algorithm_pathname, paste0("pairs_PCA_w_k", n_clusters_tanged_points_PCA, ".svg")),
        width = 8, height = 6)
    pairs(chosen_TP_w_cluster_PCA$x[,c(1:5)],
          col = tanged_points_PCA_colors[chosen_TP_w_cluster_PCA$fac$cluster],
          cex = 1,
          pch = 15,
          lower.panel = NULL)
    dev.off()
    
    
    svg(file=file.path(TP_algorithm_pathname, paste0("PCA_w_k", n_clusters_tanged_points_PCA, ".svg")),
        width = 8, height = 6)
    plot(chosen_TP_w_cluster_PCA$x[,1],
         chosen_TP_w_cluster_PCA$x[,2],
         xlab=paste("PCA 1 (", round(summary(chosen_TP_w_cluster_PCA)$importance[2]*100, 0), "%)", sep = ""),
         ylab=paste("PCA 2 (", round(summary(chosen_TP_w_cluster_PCA)$importance[5]*100, 0), "%)", sep = ""),
         bg=tanged_points_PCA_colors[chosen_TP_w_cluster_PCA$fac$cluster], 
         pch = 21,
         cex=1.5, cex.axis = 1.5, cex.lab = 1.5,
         las=1,
         asp=1)
    # Add grid lines
    abline(v=0, lty=2, col="grey50")
    abline(h=0, lty=2, col="grey50")
    dev.off()
    
    # panels for each cluster
    for (i in 1:n_clusters_tanged_points_PCA){
      
      mypath_png <- file.path(TP_algorithm_pathname_panel, paste0("chosen_TP_w_cluster_colors_cluster_", i, ".png"))
      
      png(file=mypath_png,
          width = 800, height = 800, units = "px")
      
      Momocs::panel(Momocs::slice(chosen_TP_w_cluster, cluster == i),
                    main = NULL,
                    col = tanged_points_PCA_colors[i])
      
      dev.off()
      
      
      mypath_svg <- file.path(TP_algorithm_pathname_panel, paste0("chosen_TP_w_cluster_colors_cluster_", i, ".svg"))
      
      svg(file=mypath_svg, width = 8, height = 8)
      
      Momocs::panel(Momocs::slice(chosen_TP_w_cluster, cluster == i),
                    main = NULL,
                    col = tanged_points_PCA_colors[i])
      
      dev.off()
    }
    
    
    
    # mean shapes
    min_no_of_coordinates <- list()
    for (cluster_index in 1:n_clusters_tanged_points_PCA){
      
      current_shapes <- Momocs::slice(chosen_TP_w_cluster, cluster == cluster_index)
      
      min_no_of_coordinates[[cluster_index]] <- rep(NA, length(current_shapes))
      
      for (i in 1:length(current_shapes)){
        min_no_of_coordinates[[cluster_index]][i] <- nrow(current_shapes$coo[[i]])
      }
      
      min_no_of_coordinates[[cluster_index]] <- min(min_no_of_coordinates[[cluster_index]])
      
    }
    
    ## shapes get INTERPOLATED to common number of landmarks (lowest number of landmarks per cluster)
    mean_shapes_cluster <- list()
    for (cluster_index in 1:n_clusters_tanged_points_PCA){
      mean_shapes_cluster[[cluster_index]] <- Momocs::MSHAPES(Momocs::coo_interpolate(Momocs::slice(chosen_TP_w_cluster, cluster == cluster_index),
                                                                                      n = min_no_of_coordinates[[cluster_index]])$coo)
    }
    
    chosen_TP_w_cluster_PCA_mean_shapes_cluster_out <- Momocs::Out(mean_shapes_cluster,
                                                                   fac = data.frame(cluster = paste0("cluster_", c(1:n_clusters_tanged_points_PCA))))
    
    
    
    for (i in 1:n_clusters_tanged_points_PCA){
      
      mypath_png <- file.path(TP_algorithm_pathname_means, paste0("Cluster ", i, " (n=",table(current_treecut$cluster)[[i]], ")_mean_shp.png"))
      
      png(file=mypath_png,
          width = 800, height = 800, units = "px")
      
      Momocs::panel(Momocs::slice(chosen_TP_w_cluster_PCA_mean_shapes_cluster_out, cluster == paste0("cluster_",i)),
                    main = NULL,
                    col = tanged_points_PCA_colors[i])
      
      dev.off()
      
      mypath_svg <- file.path(TP_algorithm_pathname_means, paste0("Cluster ", i, " (n=",table(current_treecut$cluster)[[i]], ")_mean_shp.svg"))
      
      svg(file=mypath_svg, width = 8, height = 8)
      
      Momocs::panel(Momocs::slice(chosen_TP_w_cluster_PCA_mean_shapes_cluster_out, cluster == paste0("cluster_",i)),
                    main = NULL,
                    col = tanged_points_PCA_colors[i])
      
      dev.off()
    }
    
    if (current_algorithm_name != "UPGMC"){
      
      # pruned tree with mean shapes
      current_algorithms_clustering_cuttree <- maptree::clip.clust(current_algorithms_clustering, 
                                                                   current_treecut, 
                                                                   k=n_clusters_tanged_points_PCA)
      ######################## change tip labels of chosen trees
      new_labels <- rep(NA, length(current_algorithms_clustering_cuttree$labels))
      for (i in 1:length(current_algorithms_clustering_cuttree$labels)){
        new_labels[i] <- paste0("Cluster ", current_algorithms_clustering_cuttree$labels[i], " (n=",table(current_treecut$cluster)[[i]], ")")
      }
    
      current_algorithms_clustering_cuttree$labels <- new_labels
      
      # 
      current_pruned_tree <- ggtree(current_algorithms_clustering_cuttree, 
                                    layout = "rectangular") + 
        xlim(NA,6) +
        geom_tiplab(aes(image=file.path(TP_algorithm_pathname_means, paste0(label, "_mean_shp.png"))), 
                    geom="image", offset=0.5, align=2, hjust = 0.5, size = 0.08) + 
        geom_tiplab(geom='label', offset=2, hjust=.5, 
                    size = 8, 
                    fontface='italic', 
                    family="TT Times New Roman") + 
        geom_treescale()
      
      ggsave(filename = file.path(TP_algorithm_pathname, paste0(current_algorithm_name, "_pruned_dendro.svg")),
             plot = current_pruned_tree,
             device = "svg",
             width = 25,
             height = 35,
             units = "cm")
    }
    
    
    
    ### distribution map
    data_wo_outliers_unique_sites <- dplyr::distinct(current_treecut, Site, .keep_all =T)
    data_wo_outliers <- current_treecut
    
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
      geom_point(data = data_wo_outliers_unique_sites[,c("Site", "Latitude", "Longitude")],  
                 aes(x = Longitude, y = Latitude, alpha = 0.9), 
                 shape = 3) +
      geom_jitter(data = data_wo_outliers,
                  aes(x = Longitude, y = Latitude,
                      fill = cluster), shape = 21, size = 3,
                  width = 0.1, height = 0.1
      ) +   
      scale_fill_manual(values = tanged_points_PCA_colors) +
      facet_wrap(~cluster,
                 scales = "fixed", 
                 labeller = as_labeller(cluster_names),
                 ncol = 3) +
      coord_quickmap() +  
      theme_classic() +  
      xlab("Longitude") +
      ylab("Latitude") +
      scale_y_continuous(expand = c(0,0)) + 
      scale_x_continuous(expand = c(0,0)) +
      theme(legend.position = "none",
            text = element_text(size=20))
    
    ggsave(filename = file.path(TP_algorithm_pathname, paste0("distribution_map_", current_algorithm_name, ".svg")),
           plot = current_distribution_map,
           device = "svg",
           width = 30,
           height = 25,
           units = "cm")
    ggsave(filename = file.path(TP_algorithm_pathname, paste0("distribution_map_", current_algorithm_name, ".png")),
           plot = current_distribution_map,
           device = "png",
           width = 30,
           height = 25,
           units = "cm")
  }
  
  
  
  
  
  ######################### NJ ##################
  
  chronozones_TP <- as.data.frame(read_csv(file.path(".","1_data","chronozones_TP.csv"), 
                                           col_types = cols(Chronozone_order = col_factor(levels = c("1", "2", "3")))))
  chronozones_TP$Chronozone <- factor(chronozones_TP$Chronozone)
  tanged_points_PCA_chronozones <- tanged_points_PCA
  tanged_points_PCA_chronozones$fac <- dplyr::left_join(tanged_points_PCA_chronozones$fac, chronozones_TP, by = "NAT")
  
  TP_NATs_NJ <- as.data.frame(tanged_points_PCA_chronozones$fac[,c("FN", "NAT", "Chronozone")])
  TP_NATs_NJ$NAT <- unlist(strsplit(TP_NATs_NJ$NAT, " Point"))
  
  chronozones_TP_df <- data.frame(tanged_points_PCA_chronozones$fac$NAT)
  rownames(chronozones_TP_df) <- tanged_points_PCA_chronozones$fac$FN
  names(chronozones_TP_df) <- "NAT"
  
  
  f <- function(x) phangorn::NJ(dist(x))
  a <- tanged_points_PCA$x[,1:minimum_no_of_pcs_tanged_points_PCA]
  tr <- f(a)
  X <- ape::boot.phylo(phy = tr, 
                       x = a, 
                       FUN = f, 
                       trees = TRUE,
                       B = 10000,
                       multicore = T,
                       mc.cores = parallel::detectCores()-1) 
  saveRDS(X,
          file = file.path(output_folder, "NJ_bootstrapped10k_trees_TP.RDS"))
  
  tree <- phangorn::plotBS(tr, X$trees)
  saveRDS(tree,
          file = file.path(output_folder, "NJ_bootstrapped10k_trees_TP_tree.RDS"))
  tree <- readRDS(file.path(output_folder, "NJ_bootstrapped10k_trees_TP_tree.RDS"))
  tree2 <- phangorn::pruneTree(tree, 50)#mayority consensus tree
  
  
  tanged_points_PCA_NJ <- tree
  tanged_points_PCA_NJ$tip.label[which(tanged_points_PCA_NJ$tip.label=="H\\E4cklingen_11")] <- "Haecklingen_11"
  
  NJ_TP_ggtree <- 
    ggtree(tanged_points_PCA_NJ) %<+% TP_NATs_NJ +
    geom_tiplab(size=5, 
                aes(color = Chronozone),
                align = F) +
    geom_treescale() + 
    scale_colour_discrete(na.translate = F) #+ 
    # # guides(fill=guide_legend(title="NAT")) +
    # geom_label2(aes(subset = (!isTip & as.numeric(label) >= 50),
    #                 label = round(as.numeric(label), digits = 0))) # +
  # geom_label2(aes(subset = (!isTip), label=round(as.numeric(label), digits = 0)))
  
  
  
  chronozones_TP_heatmap <- data.frame(NAT = TP_NATs_NJ$NAT)
  rownames(chronozones_TP_heatmap) <- TP_NATs_NJ$FN
  
  NJ_TP_gheatmap <- gheatmap(NJ_TP_ggtree,
                             chronozones_TP_heatmap,
                             offset=0.1,
                             width=0.1,
                             colnames = T,
                             font.size = 10) +
    scale_x_ggtree() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.text=element_text(size=31),
          legend.title=element_text(size=31)) + 
    guides(fill=guide_legend(title="NAT")) +
    scale_fill_discrete(na.translate = F)
  
    NJ_TP_gheatmap_height <- 30
    NJ_TP_gheatmap_width <- 15
    
    ggsave(NJ_TP_gheatmap,
           filename = file.path(output_folder, "Fig_10_NJ_bootstrapped10k_with_heatmap_chronozones.svg"),
           width = NJ_TP_gheatmap_width,
           height = NJ_TP_gheatmap_height,
           units = "in")
    ggsave(NJ_TP_gheatmap,
           filename = file.path(output_folder, "Fig_10_NJ_bootstrapped10k_with_heatmap_chronozones.png"),
           width = NJ_TP_gheatmap_width,
           height = NJ_TP_gheatmap_height,
           units = "in")
      
  