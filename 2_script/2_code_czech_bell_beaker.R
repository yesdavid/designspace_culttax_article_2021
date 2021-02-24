library(cluster)
library(ggplot2)
library(readr)
library(Momocs)
library(dplyr)

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

output_folder <- file.path(".","3_output","czech_bell_beaker_arrowheads")
dir.create(output_folder,
           recursive = T)

# load outlines
outlines_combined_petrik <- readRDS(file = file.path(".","1_data","outlines_combined_petrik_2018.RDS"))


# GMM procedures
outlines_combined_petrik_centered <- Momocs::coo_centre(outlines_combined_petrik) # center
outlines_combined_petrik_centered_scaled <- Momocs::coo_scale(outlines_combined_petrik_centered) # scale
outlines_combined_petrik_centered_scaled <- Momocs::coo_slidedirection(outlines_combined_petrik_centered_scaled,
                                                                       direction = "up") 
# outline inspection
# stack(outlines_combined_petrik_centered_scaled)
# Momocs::panel(outlines_combined_petrik_centered_scaled,
#               col = "grey")

# harmonic calibration. Estimates the number of harmonics required for the Fourier methods implemented in Momocs. This is the only step in this section that produces data we need in the subsequent step.
outlines_combined_petrik_centered_scaled_harmonics <- Momocs::calibrate_harmonicpower_efourier(outlines_combined_petrik_centered_scaled, 
                                                                                               plot = F)
outlines_combined_petrik_centered_scaled_harmonics

# efourier
outlines_combined_petrik_centered_scaled_efourier <- Momocs::efourier(outlines_combined_petrik_centered_scaled,
                                                                      nb.h = as.matrix(outlines_combined_petrik_centered_scaled_harmonics[["minh"]])[[4,1]], # choses number of harmonics for 99.9%
                                                                      norm = F) 
# PCA
outlines_combined_petrik_centered_scaled_PCA <- Momocs::PCA(outlines_combined_petrik_centered_scaled_efourier) # PCA on Coe objects, using prcomp.

minimum_no_of_pcs_petrik <- Momocs::scree_min(outlines_combined_petrik_centered_scaled_PCA,
                                              prop = 0.95) 
minimum_no_of_pcs_petrik

petrik_screeplot <- Momocs::scree_plot(outlines_combined_petrik_centered_scaled_PCA,
                                       nax = 1:8)
ggsave(petrik_screeplot,
       filename = file.path(output_folder, "petrik_screeplot.png"),
       width = 5,
       height = 5)

# PCA shape variation
gg <- Momocs::PCcontrib(outlines_combined_petrik_centered_scaled_PCA,
                        nax = 1:minimum_no_of_pcs_petrik,
                        sd.r = c(-2,-1,0,1,2)) 

petrik_pccontrib <- gg$gg + 
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

ggsave(petrik_pccontrib,
       filename = file.path(output_folder, "petrik_pccontrib.png"),
       width = 5,
       height = 5)

# PCA plot
plot(outlines_combined_petrik_centered_scaled_PCA, 
     xax = 1, 
     yax = 2, 
     pos.shp = "XY",
     rug = FALSE, 
     zoom = 1, 
     lwd.shp = 1, 
     size.shp = 0.5, 
     amp.shp = 0.5, 
     title = "")


###

pairs(outlines_combined_petrik_centered_scaled_PCA$x[,c(1:5)],
      lower.panel = NULL)

###
# UPGMA (unweighted pair-group method using arithmetic average)
outlines_combined_petrik_centered_scaled_PCA_as_matrix <- as.matrix(outlines_combined_petrik_centered_scaled_PCA$x[,1:minimum_no_of_pcs_petrik])
petrik_dist <- dist(outlines_combined_petrik_centered_scaled_PCA_as_matrix)

outlines_petrik_upgma <- hclust(petrik_dist,
                                method = "average")

library(NbClust)
petrik_NbClust_ward <- NbClust::NbClust(data = outlines_combined_petrik_centered_scaled_PCA$x[,1:minimum_no_of_pcs_petrik],
                                        distance = "euclidean",
                                        method = "average",
                                        index = c("gap", "silhouette"))
petrik_NbClust <- petrik_NbClust_ward$All.index
petrik_NbClust_df <- as.data.frame(petrik_NbClust)
petrik_NbClust_df$NClust <- 1:nrow(petrik_NbClust)+1
silhouette_plot <- ggplot(petrik_NbClust_df, aes(x = NClust, y = Silhouette)) + geom_point() + geom_line() + theme_bw() +
  scale_x_continuous(breaks = seq(2, max(petrik_NbClust_df$NClust)+1, by = 1)) +
  xlab("Number of clusters") +
  ylab("Average silhouette value")
silhouette_plot
gap_plot <- ggplot(petrik_NbClust_df, aes(x = NClust, y = Gap)) + geom_point() + geom_line() + theme_bw() +
  scale_x_continuous(breaks = seq(2, max(petrik_NbClust_df$NClust)+1, by = 1)) +
  xlab("Number of clusters") +
  ylab("Gap statistic (k)")
gap_plot
# cowplot::plot_grid(silhouette_plot, gap_plot)

height_silhouette_plot <- 100
width_silhouette_plot <- 100

ggsave(silhouette_plot,
       filename = file.path(output_folder, "silhouette_plot_NbClust_average.svg"),
       height  = height_silhouette_plot,
       width = width_silhouette_plot,
       units = "mm")
ggsave(silhouette_plot,
       filename = file.path(output_folder, "silhouette_plot_NbClust_average.png"),
       height  = height_silhouette_plot,
       width = width_silhouette_plot,
       units = "mm")


petrik_colors <- c("darkgreen","darkblue","maroon3","red","yellow","burlywood","chartreuse","cyan","deeppink","deepskyblue1")
petrik_colors_code <- c("#006400","#00008b","#b03060","#ff0000","#ffff00","#deb887","#7FFF00","#00ffff","#ff00ff","#6495ed")


n_clusters_petrik <- 10


###


clusters_petrik_cutree_k10 <- data.frame(ID_artefact = row.names(as.data.frame(cutree(outlines_petrik_upgma, k = n_clusters_petrik))),
                                         cluster = as.factor(as.data.frame(cutree(outlines_petrik_upgma, k = n_clusters_petrik))[[1]]))

outlines_combined_petrik_centered_scaled_w_cluster <- Momocs::Out(outlines_combined_petrik_centered_scaled$coo,
                                                                  fac = clusters_petrik_cutree_k10)



# # saves plot of each individual cluster with its respective artefacts
upgma_path <- file.path(output_folder "petrik_clusters_panels_upgma")
dir.create(upgma_path,
           recursive = T)

for (i in 1:n_clusters_petrik){

  mypath <- file.path(upgma_path, paste0("petrik_2018_upgma_k10_cluster_colors_cluster_", i, ".png"))

  png(file=mypath,
      width = 800, height = 800, units = "px")
  Momocs::panel(Momocs::slice(outlines_combined_petrik_centered_scaled_w_cluster, cluster == i),
                # main = paste("Cluster", i, "; n = ", nrow(filter(outlines_combined_petrik_centered_scaled_w_cluster$fac, cluster == i))),
                col = petrik_colors[i])
  dev.off()
}


# number of artefacts per cluster
petrik_2016_upgma_k10_cluster_no_of_artefacts <- list()
for (i in 1:n_clusters_petrik){
  petrik_2016_upgma_k10_cluster_no_of_artefacts[[i]] <- nrow(subset(outlines_combined_petrik_centered_scaled_w_cluster$fac, cluster == i))
}

number_of_artefacts_per_cluster_petrik <- data.frame(cluster = c(1:n_clusters_petrik),
                                                     n = do.call(rbind, petrik_2016_upgma_k10_cluster_no_of_artefacts)[,1]) # number of artefacts per cluster
nrow(outlines_combined_petrik_centered_scaled_w_cluster$fac) # total number of artefacts


outlines_combined_petrik_centered_scaled_w_cluster_harmonics <- Momocs::calibrate_harmonicpower_efourier(outlines_combined_petrik_centered_scaled_w_cluster,plot = F)
outlines_combined_petrik_centered_scaled_w_cluster_efourier <- Momocs::efourier(outlines_combined_petrik_centered_scaled_w_cluster,
                                                                                nb.h = as.matrix(outlines_combined_petrik_centered_scaled_w_cluster_harmonics[["minh"]])[[4,1]], norm = F) 
outlines_combined_petrik_centered_scaled_w_cluster_PCA <- Momocs::PCA(outlines_combined_petrik_centered_scaled_w_cluster_efourier) # PCA on Coe objects, using prcomp.



## PCA plot ggplot
petrik_w_Cluster_PCA_names <- as.data.frame(outlines_combined_petrik_centered_scaled_w_cluster_PCA$x[,c(1:minimum_no_of_pcs_petrik)])
rownames(petrik_w_Cluster_PCA_names) <- outlines_combined_petrik_centered_scaled_w_cluster_PCA$fac$ID_artefact
petrik_w_Cluster_PCA_names$Cluster <- outlines_combined_petrik_centered_scaled_w_cluster_PCA$fac$cluster

names(petrik_colors_code) <- c(1:length(unique(outlines_combined_petrik_centered_scaled_w_cluster_PCA$fac$cluster)))

petrik_w_Cluster_PCA_names$color <- sapply(petrik_w_Cluster_PCA_names$Cluster, function(x){petrik_colors_code[as.integer(x)]})

PCA_plot_petrik_with_outliers <- ggplot(data = petrik_w_Cluster_PCA_names, aes(x = PC1, y = PC2, fill = Cluster)) +
  geom_point(size = 3, pch = 21) +
  scale_fill_manual(values = petrik_colors_code) +
  coord_fixed(ratio =1) +
  theme_bw() +
  xlab(paste0("PC1 (", round(outlines_combined_petrik_centered_scaled_w_cluster_PCA$eig[1]*100, digits = 0), "%)")) +
  ylab(paste0("PC2 (", round(outlines_combined_petrik_centered_scaled_w_cluster_PCA$eig[2]*100, digits = 0), "%)")) +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)) + 
  guides(fill = FALSE) + 
  geom_hline(yintercept=0, linetype="dashed", alpha = 0.5) + 
  geom_vline(xintercept=0, linetype="dashed", alpha = 0.5)


ggsave(PCA_plot_petrik_with_outliers,
       filename = file.path(output_folder, "PCA_plot_petrik_with_outliers.svg"),
       width = 7,
       height = 5,
       dpi = 320,
       units = "in",
       device = "svg")


# mean shapes
min_no_of_coordinates <- list()
for (cluster_index in 1:n_clusters_petrik){
  current_shapes <- Momocs::slice(outlines_combined_petrik_centered_scaled_w_cluster, cluster == cluster_index)
  min_no_of_coordinates[[cluster_index]] <- rep(NA, length(current_shapes))
  for (i in 1:length(current_shapes)){
    min_no_of_coordinates[[cluster_index]][i] <- nrow(current_shapes$coo[[i]])
  }
  min_no_of_coordinates[[cluster_index]] <- min(min_no_of_coordinates[[cluster_index]])
}

# shapes get INTERPOLATED to common number of landmarks (lowest number of landmarks per cluster)
mean_shapes_petrik_cluster <- list()
for (cluster_index in 1:n_clusters_petrik){
  mean_shapes_petrik_cluster[[cluster_index]] <- Momocs::MSHAPES(Momocs::coo_interpolate(Momocs::slice(outlines_combined_petrik_centered_scaled_w_cluster, cluster == cluster_index),
                                                                                         n = min_no_of_coordinates[[cluster_index]])$coo)
}


mean_shapes_petrik_cluster_out <- Momocs::Out(mean_shapes_petrik_cluster,
                                              fac = data.frame(cluster = paste0("cluster_", c(1:n_clusters_petrik))))

panel(mean_shapes_petrik_cluster_out, 
      names = T,
      main = "Mean shapes of each UPGMA cluster (Petrík et al. 2018 data as outlines)",
      cex.names = 1.5,
      col = "grey")


####### inspect outliers from PCA; all axis

pairs(outlines_combined_petrik_centered_scaled_w_cluster_PCA$x[,c(1:minimum_no_of_pcs_petrik)],
      lower.panel = NULL, 
      col = petrik_colors[outlines_combined_petrik_centered_scaled_w_cluster_PCA$fac$cluster], 
      cex = 1, 
      pch = 15)



############## create and plot dendrogram with mean shapes at the tips ############## 

library(cluster)
library(maptree)
library(ggimage)
library(ggtree)

outlines_petrik_upgma_cuttree <- maptree::clip.clust(outlines_petrik_upgma, 
                                                     clusters_petrik_cutree_k10, 
                                                     k=n_clusters_petrik)

# mean shapes

mean_shapes_upgma_path <- file.path(output_folder, "petrik_clusters_mean_shapes_upgma")
dir.create(mean_shapes_upgma_path,
           recursive = T)

for (i in 1:n_clusters_petrik){
  
  mypath <- file.path(mean_shapes_upgma_path, paste0("petrik_2018_w_cluster_upgma_colors_", "Cluster ", i, " (n=",colSums(table(clusters_petrik_cutree_k10))[[i]], ")_mean_shp.png"))
  
  png(file=mypath,
      width = 800, height = 800, units = "px")
  Momocs::panel(Momocs::slice(mean_shapes_petrik_cluster_out, cluster == paste0("cluster_",i)),
                main = NULL,
                col = petrik_colors[i])
  
  dev.off()
}

layout = "circular"
img_size <- 0.02
tip_lab_size <- 3
path_to_mean_shapes_petrik <- file.path(mean_shapes_upgma_path, "petrik_2018_w_cluster_upgma_colors_")




######################## change tip labels of chosen trees
new_labels <- rep(NA, length(outlines_petrik_upgma_cuttree$labels))
for (i in 1:length(outlines_petrik_upgma_cuttree$labels)){
  new_labels[i] <- paste0("Cluster ", i, " (n=",colSums(table(clusters_petrik_cutree_k10))[[i]], ")")
}

outlines_petrik_upgma_cuttree$labels <- new_labels


petrik_pruned_tree <- ggtree(outlines_petrik_upgma_cuttree, 
                             layout = "rectangular") + 
  xlim(NA, 6) + 
  geom_tiplab(aes(image=paste0(path_to_mean_shapes_petrik, label, "_mean_shp.png")), 
              geom="image", offset=2, align=2, size=0.1, hjust = 0.5)  + 
  geom_tiplab(geom='label', offset=1, hjust=.5, size = 5) + 
  geom_treescale()

ggsave(filename = file.path(output_folder, paste0("pruned_dendro_UPGMA_petrik_k_", n_clusters_petrik, ".svg")),
       plot = petrik_pruned_tree,
       device = "svg",
       width = 25,
       height = 35,
       units = "cm")











#########################################################################################################
################################### with outliers removed ###############################################
#########################################################################################################

# Method D): DBSCAN
petrik_outliers_db <- fpc::dbscan(outlines_combined_petrik_centered_scaled_PCA$x, eps = 0.3, MinPts = 3)
plot(petrik_outliers_db, outlines_combined_petrik_centered_scaled_PCA$x, main = "DBSCAN", frame = FALSE)


petrik_outliers_Cluster <- data.frame(ID_artefact = row.names(outlines_combined_petrik_centered_scaled_PCA$x), 
                                      Cluster = petrik_outliers_db$cluster,
                                      row.names = NULL)

petrik_outliers_Cluster_outlier_names <- subset(petrik_outliers_Cluster, Cluster != 1)

outlines_combined_petrik_centered_scaled_PCA$fac$outlier_names <- NA

outlines_combined_petrik_centered_scaled_w_outliersCluster <- Momocs::Out(outlines_combined_petrik_centered_scaled$coo,
                                                                          fac = petrik_outliers_Cluster)


# without outliers
petrik_without_outliers <- Momocs::slice(outlines_combined_petrik_centered_scaled, -match(petrik_outliers_Cluster_outlier_names$ID_artefact, outlines_combined_petrik_centered_scaled_w_outliersCluster$fac$ID_artefact))




# harmonic calibration. Estimates the number of harmonics required for the Fourier methods implemented in Momocs. This is the only step in this section that produces data we need in the subsequent step.
petrik_without_outliers_harmonics <- Momocs::calibrate_harmonicpower_efourier(petrik_without_outliers, 
                                                                              plot = F)

# efourier
petrik_without_outliers_efourier <- Momocs::efourier(petrik_without_outliers,
                                                     nb.h = as.matrix(petrik_without_outliers_harmonics[["minh"]])[[4,1]], # choses number of harmonics for 99.9%
                                                     norm = F) 
# PCA
petrik_without_outliers_PCA <- Momocs::PCA(petrik_without_outliers_efourier) # PCA on Coe objects, using prcomp.

minimum_no_of_pcs_petrik <- Momocs::scree_min(petrik_without_outliers_PCA,
                                              prop = 0.95) 
minimum_no_of_pcs_petrik

Momocs::scree_plot(petrik_without_outliers_PCA,
                   nax = 1:(minimum_no_of_pcs_petrik+2))

minimum_no_of_pcs_petrik <- 6

# PCA shape variation
gg <- Momocs::PCcontrib(petrik_without_outliers_PCA,
                        nax = 1:minimum_no_of_pcs_petrik, # sind alle minimum_no_of_pcs_petrik PCs "signifikant"/beschreiben alle PCs eine morphologische Eigenschaft? 
                        sd.r = c(-2,0,2)) 

# PCA plot
plot(petrik_without_outliers_PCA, 
     xax = 1, 
     yax = 2, 
     pos.shp = "XY",
     rug = FALSE, 
     zoom = 1, 
     lwd.shp = 1, 
     size.shp = 0.5, 
     amp.shp = 0.5, 
     title = "")



#### Cluster analysis

library(NbClust)
petrik_NbClust_ward <- NbClust::NbClust(data = petrik_without_outliers_PCA$x,
                                        distance = "euclidean",
                                        method = "average",
                                        index = c("gap", "silhouette"))
petrik_NbClust <- petrik_NbClust_ward$All.index
petrik_NbClust_df <- as.data.frame(petrik_NbClust)
petrik_NbClust_df$NClust <- 1:nrow(petrik_NbClust)+1
ggplot(petrik_NbClust_df, aes(x = NClust, y = Silhouette)) + geom_point() + geom_line() + theme_bw() +
  scale_x_continuous(breaks = seq(2, max(petrik_NbClust_df$NClust)+1, by = 1)) +
  xlab("Number of clusters") +
  ylab("Average silhouette value")





petrik_without_outliers_colors <- c("darkgreen","darkblue","maroon3","red","yellow","burlywood","chartreuse","cyan","deeppink","deepskyblue1")
petrik_without_outliers_colors_code <- c("#006400","#00008b","#b03060","#ff0000","#ffff00","#deb887","#7FFF00","#00ffff","#ff00ff","#6495ed")


n_Clusters_petrik_without_outliers <- 6

outlines_petrik_without_outliers_upgma <- hclust(dist(petrik_without_outliers_PCA$x),
                                                 method = "average")

###


Clusters_petrik_without_outliers_cutree <- data.frame(ID_artefact = row.names(as.data.frame(cutree(outlines_petrik_without_outliers_upgma, 
                                                                                                   k = n_Clusters_petrik_without_outliers))),
                                                      Cluster = as.factor(as.data.frame(cutree(outlines_petrik_without_outliers_upgma, k = n_Clusters_petrik_without_outliers))[[1]]))

petrik_without_outliers_w_Cluster <- Momocs::Out(petrik_without_outliers$coo,
                                                 fac = Clusters_petrik_without_outliers_cutree)





# saves plot of each individual Cluster with its respective artefacts

path_no_outliers_upgma <- file.path(output_folder, "petrik_without_outliers_clusters_panels_upgma")
dir.create(path_no_outliers_upgma,
           recursive = T)

for (i in 1:n_Clusters_petrik_without_outliers){
  
  mypath_png <- file.path(path_no_outliers_upgma, 
                       paste0("petrik_without_outliers_2018_upgma_k10_cluster_colors_cluster_", i, 
                       ".png"))
  mypath_svg <- file.path(path_no_outliers_upgma, 
                       paste0("petrik_without_outliers_2018_upgma_k10_cluster_colors_cluster_", i, 
                       ".svg"))
  
  png(file=mypath_png,
      width = 800, height = 800, units = "px")
  Momocs::panel(Momocs::slice(petrik_without_outliers_w_Cluster, Cluster == i),
                col = petrik_without_outliers_colors[i])
  dev.off()
  svg(file=mypath_svg,
      width = 8, height = 8)
  Momocs::panel(Momocs::slice(petrik_without_outliers_w_Cluster, Cluster == i),
                col = petrik_without_outliers_colors[i])
  dev.off()
}

# plots each individual Cluster with its respective artefacts
for (i in 1:n_Clusters_petrik_without_outliers){
  Momocs::panel(Momocs::slice(petrik_without_outliers_w_Cluster, Cluster == i),
                main = paste("Cluster", i, "; n = ", nrow(subset(petrik_without_outliers_w_Cluster$fac, Cluster == i))),
                col = petrik_without_outliers_colors[i])
}

# number of artefacts per Cluster
petrik_without_outliers_2016_upgma_k10_Cluster_no_of_artefacts <- list()
for (i in 1:n_Clusters_petrik_without_outliers){
  petrik_without_outliers_2016_upgma_k10_Cluster_no_of_artefacts[[i]] <- nrow(subset(petrik_without_outliers_w_Cluster$fac, Cluster == i))
}

number_of_artefacts_per_Cluster_petrik_without_outliers <- data.frame(Cluster = c(1:n_Clusters_petrik_without_outliers),
                                                                      n = do.call(rbind, petrik_without_outliers_2016_upgma_k10_Cluster_no_of_artefacts)[,1]) # number of artefacts per Cluster
nrow(petrik_without_outliers_w_Cluster$fac) # total number of artefacts

petrik_without_outliers_w_Cluster_harmonics <- Momocs::calibrate_harmonicpower_efourier(petrik_without_outliers_w_Cluster,plot = F)
petrik_without_outliers_w_Cluster_efourier <- Momocs::efourier(petrik_without_outliers_w_Cluster,
                                                               nb.h = as.matrix(petrik_without_outliers_w_Cluster_harmonics[["minh"]])[[4,1]], norm = F) 
petrik_without_outliers_w_Cluster_PCA <- Momocs::PCA(petrik_without_outliers_w_Cluster_efourier) # PCA on Coe objects, using prcomp.



## PCA plot ggplot
petrik_without_outliers_w_Cluster_PCA_names <- as.data.frame(petrik_without_outliers_w_Cluster_PCA$x[,c(1:minimum_no_of_pcs_petrik)])
rownames(petrik_without_outliers_w_Cluster_PCA_names) <- petrik_without_outliers_w_Cluster_PCA$fac$ID_artefact
petrik_without_outliers_w_Cluster_PCA_names$Cluster <- petrik_without_outliers_w_Cluster_PCA$fac$Cluster

names(petrik_without_outliers_colors_code) <- c(1:length(unique(petrik_without_outliers_w_Cluster_PCA$fac$Cluster)))

petrik_without_outliers_w_Cluster_PCA_names$color <- sapply(petrik_without_outliers_w_Cluster_PCA_names$Cluster, function(x){petrik_without_outliers_colors_code[as.integer(x)]})

PCA_plot_petrikWOoutliers <- ggplot(data = petrik_without_outliers_w_Cluster_PCA_names, aes(x = PC1, y = PC2, fill = Cluster)) +
  geom_point(size = 3, pch = 21) +
  scale_fill_manual(values = petrik_without_outliers_colors_code) +
  coord_fixed(ratio =1) +
  theme_bw() +
  xlab(paste0("PC1 (", round(petrik_without_outliers_w_Cluster_PCA$eig[1]*100, digits = 0), "%)")) +
  ylab(paste0("PC2 (", round(petrik_without_outliers_w_Cluster_PCA$eig[2]*100, digits = 0), "%)")) +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)) + 
  guides(fill = FALSE) + 
  geom_hline(yintercept=0, linetype="dashed", alpha = 0.5) + 
  geom_vline(xintercept=0, linetype="dashed", alpha = 0.5)


ggsave(PCA_plot_petrikWOoutliers,
       filename = file.path(output_folder, "PCA_plot_petrikWOoutliers.svg"),
       width = 7,
       height = 5,
       dpi = 320,
       units = "in",
       device = "svg")

# mean shapes
min_no_of_coordinates <- list()
for (Cluster_index in 1:n_Clusters_petrik_without_outliers){
  current_shapes <- Momocs::slice(petrik_without_outliers_w_Cluster, Cluster == Cluster_index)
  min_no_of_coordinates[[Cluster_index]] <- rep(NA, length(current_shapes))
  for (i in 1:length(current_shapes)){
    min_no_of_coordinates[[Cluster_index]][i] <- nrow(current_shapes$coo[[i]])
  }
  min_no_of_coordinates[[Cluster_index]] <- min(min_no_of_coordinates[[Cluster_index]])
}

# shapes get INTERPOLATED to common number of landmarks (lowest number of landmarks per Cluster)
mean_shapes_petrik_without_outliers_Cluster <- list()
for (Cluster_index in 1:n_Clusters_petrik_without_outliers){
  mean_shapes_petrik_without_outliers_Cluster[[Cluster_index]] <- Momocs::MSHAPES(Momocs::coo_interpolate(Momocs::slice(petrik_without_outliers_w_Cluster, Cluster == Cluster_index),
                                                                                                          n = min_no_of_coordinates[[Cluster_index]])$coo)
}


mean_shapes_petrik_without_outliers_Cluster_out <- Momocs::Out(mean_shapes_petrik_without_outliers_Cluster,
                                                               fac = data.frame(Cluster = paste0("Cluster_", c(1:n_Clusters_petrik_without_outliers))))

panel(mean_shapes_petrik_without_outliers_Cluster_out, 
      names = T,
      main = "Mean shapes of each UPGMA Cluster (Petrík et al. 2018 data as outlines)",
      cex.names = 1.5,
      col = "grey")




outlines_petrik_without_outliers_upgma_cuttree <- maptree::clip.clust(outlines_petrik_without_outliers_upgma, 
                                                                      Clusters_petrik_without_outliers_cutree, 
                                                                      k=n_Clusters_petrik_without_outliers)

# mean shapes
path_petrik_without_outliers_clusters_mean_shapes_upgma <- file.path(output_folder, "petrik_without_outliers_clusters_mean_shapes_upgma")
dir.create(path_petrik_without_outliers_clusters_mean_shapes_upgma,
           recursive = T)

for (i in 1:n_Clusters_petrik_without_outliers){
  
  mypath_png <- file.path(path_petrik_without_outliers_clusters_mean_shapes_upgma,
                          paste0("petrik_without_outliers_2018_w_cluster_upgma_colors_Cluster ", 
                              i, 
                              " (n=",
                              colSums(table(Clusters_petrik_without_outliers_cutree))[[i]], 
                              ")_mean_shp.png"))
  mypath_svg <- file.path(path_petrik_without_outliers_clusters_mean_shapes_upgma,
                          paste0("petrik_without_outliers_2018_w_cluster_upgma_colors_Cluster ", 
                                 i, 
                                 " (n=",
                                 colSums(table(Clusters_petrik_without_outliers_cutree))[[i]], 
                                 ")_mean_shp.svg"))
  
  png(file=mypath_png,
      width = 800, height = 800, units = "px")
  Momocs::panel(Momocs::slice(mean_shapes_petrik_without_outliers_Cluster_out, Cluster == paste0("Cluster_",i)),
                main = NULL,
                col = petrik_without_outliers_colors[i])
  dev.off()
  
  svg(file=mypath_svg,
      width = 8, height = 8)
  Momocs::panel(Momocs::slice(mean_shapes_petrik_without_outliers_Cluster_out, Cluster == paste0("Cluster_",i)),
                main = NULL,
                col = petrik_without_outliers_colors[i])
  dev.off()
}

layout = "circular"
img_size <- 0.02
tip_lab_size <- 3
path_to_mean_shapes_petrik_without_outliers <- file.path(path_petrik_without_outliers_clusters_mean_shapes_upgma, "petrik_without_outliers_2018_w_cluster_upgma_colors_")



######################## change tip labels of chosen trees
new_labels <- rep(NA, length(outlines_petrik_without_outliers_upgma_cuttree$labels))
for (i in 1:length(outlines_petrik_without_outliers_upgma_cuttree$labels)){
  new_labels[i] <- paste0("Cluster ", i, " (n=",colSums(table(Clusters_petrik_without_outliers_cutree))[[i]], ")")
}

outlines_petrik_without_outliers_upgma_cuttree$labels <- new_labels


petrik_without_outliers_pruned_tree <- ggtree(outlines_petrik_without_outliers_upgma_cuttree, 
                                              layout = "rectangular") + 
  xlim(NA, 6) + 
  geom_tiplab(aes(image=paste0(path_to_mean_shapes_petrik_without_outliers, label, "_mean_shp.png")), 
              geom="image", offset=0.1, align=2, size=0.1, hjust = 0.5)  + 
  geom_tiplab(geom='label', offset=0.05, hjust=.5, size = 5) + 
  geom_treescale() +
  xlim(0,0.2)

ggsave(filename = file.path(output_folder,
                            paste0("pruned_dendro_UPGMA_petrik_without_outliers_k_", 
                                   n_Clusters_petrik_without_outliers, 
                                   ".svg")),
       plot = petrik_without_outliers_pruned_tree,
       device = "svg",
       width = 25,
       height = 35,
       units = "cm")



