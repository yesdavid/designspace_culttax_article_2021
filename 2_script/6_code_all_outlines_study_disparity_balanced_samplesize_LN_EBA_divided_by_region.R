# overall variability/disparity

library(ggplot2)
library(readr)
library(fpc)
library(Momocs)
library(dplyr)

rm(list=ls())


#### BUCHANAN ET AL.
# load outlines
outlines_combined_buchanan <- readRDS(file.path(".", "1_data", "outlines_combined_goshen_plainview.RDS"))
# rename artefact names to match spreadsheet
buchanan_names <- list()
for (i in 1:length(outlines_combined_buchanan$coo)){
  buchanan_names[i] <- strsplit(names(outlines_combined_buchanan$coo), split = "rotated-")[[i]][2]
}
names(outlines_combined_buchanan$coo) <- as.vector(do.call(rbind, buchanan_names))
# load spreadsheet
buchanan_supp <- as.data.frame(readr::read_csv(file.path(".", "1_data", "Buchanan_et_al_AA_Supp_Mats_edited.csv"),
                                               col_types = cols(prev_assigned_type = col_factor(levels = c("Goshen", "Plainview")))))
# add spreadsheet data to outlines
outlines_combined_buchanan <- Momocs::Out(outlines_combined_buchanan$coo,
                                          fac = buchanan_supp)


#### PETRIK ET AL.
# load outlines
outlines_combined_petrik <- readRDS(file = file.path(".", "1_data", "outlines_combined_petrik_2018.RDS"))


#### NICOLAS
# load outlines
outlines_combined_nicolas <- readRDS(file = file.path(".", "1_data", "outlines_combined_nicholas_2016.RDS"))

# unification of outlines with catalogue-dataframe
nicolas_fleches_2016_catalog_ids_coordinates <- readr::read_csv(file = file.path(".", "1_data", "nicolas_fleches_2016_catalog_ids_with_coordinates.csv"))

outlines_combined_nicolas_2016_names <- names(outlines_combined_nicolas)
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



outlines_combined_nicolas_fac <- Momocs::Out(outlines_combined_nicolas$coo,
                                             fac = nicolas_fleches_2016_catalog_ids_coordinates_artefact_ID)


outlines_nicolas_dk <- Momocs::filter(outlines_combined_nicolas_fac,
                                      country %in% "Denmark")
outlines_nicolas_uk <- Momocs::filter(outlines_combined_nicolas_fac,
                                      country %in% "United Kingdom")
outlines_nicolas_fr <- Momocs::filter(outlines_combined_nicolas_fac,
                                      country %in% "France")


#### Tanged points
# load outlines
tanged_points_tps <- Momocs::import_tps(file.path(".", "1_data", "TPS_TP_27_09_2019.TPS"), 
                                        curves = TRUE)

tanged_points <- readr::read_csv(file.path(".", "1_data", "tanged.points.csv"),
                                 col_types = cols(Context = col_character(),
                                                  Site = col_character()))
# subset
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
outlinetpNAT <- Momocs::filter(outlinetp, 
                               !NAT %in% c("Tanged Point (Unspecified)", NA))
outlinetpNAT_subset  <- Momocs::filter(outlinetpNAT, 
                                       !NAT %in% c("Swiderian Point", "Chwalibogowice Point", "Hamburgian Point")) # we remove those NAT-categories from the dataset.

table(outlinetpNAT_subset$fac$Country)
table(outlinetpNAT_subset$fac$NAT)






################# COMBINE OUTLINES

names_buchanan_df <- data.frame(ID = names(outlines_combined_buchanan),
                                Dataset = "Palaeoindian")
names_petrik_df <- data.frame(ID = names(outlines_combined_petrik),
                              Dataset = "Bell Beaker CZ")
names_nicolas_fr <- data.frame(ID = names(outlines_nicolas_fr),
                               Dataset = "Late Neolithic/\nEarly Bronze Age\nFR")
names_nicolas_uk <- data.frame(ID = names(outlines_nicolas_uk),
                               Dataset = "Late Neolithic/\nEarly Bronze Age\nUK")
names_nicolas_dk <- data.frame(ID = names(outlines_nicolas_dk),
                               Dataset = "Late Neolithic/\nEarly Bronze Age\nDK")
names_TP_df <- data.frame(ID = names(outlinetpNAT_subset),
                          Dataset = "Tanged points")

outlines_bellbeaker_all <- rbind(names_petrik_df, names_nicolas_fr, names_nicolas_uk, names_nicolas_dk)

all_outlines_df <- rbind(names_buchanan_df, names_petrik_df, names_nicolas_fr, names_nicolas_uk, names_nicolas_dk, names_TP_df) 

all_outlines <- c(outlines_combined_buchanan$coo, outlines_combined_petrik$coo, outlines_nicolas_fr$coo,  outlines_nicolas_uk$coo, outlines_nicolas_dk$coo, outlinetpNAT_subset$coo) 

row.names(all_outlines_df) <- all_outlines_df$ID
all_outlines_OUT <- Momocs::Out(all_outlines,
                                fac = all_outlines_df)


##### try to balance out the distribution between the three datasets. to fit the TP dataset in size

# library(splitstackshape)
# set.seed(1)
# stratified_sample_BellBeaker <- splitstackshape::stratified(outlines_bellbeaker_all, 
#                                                  group = "Dataset", 
#                                                  size = nrow(names_TP_df)/length(unique(outlines_bellbeaker_all$Dataset)), 
#                                                  replace = F,
#                                                  bothSets = F)
# 
# all_outlines_df_resampled <- as.data.frame(rbind(names_buchanan_df, stratified_sample_BellBeaker, names_TP_df))



## in this case (Late Neolithic/Early Bronze Age), the data gets resampled to the lowest number of artefacts available in table(outlines_bellbeaker_all$Dataset)
library(splitstackshape)
set.seed(1)
stratified_sample_BellBeaker <- splitstackshape::stratified(outlines_bellbeaker_all,
                                                 group = "Dataset",
                                                 size = min(table(outlines_bellbeaker_all$Dataset)),
                                                 replace = F,
                                                 bothSets = F)

all_outlines_df_resampled <- as.data.frame(rbind(names_buchanan_df, stratified_sample_BellBeaker, names_TP_df))


resampled_shapes <- Momocs::filter(all_outlines_OUT,
                                   ID %in% all_outlines_df_resampled$ID)

# resampled_shapes$fac$Dataset <- as.character(resampled_shapes$fac$Dataset)


all_outlines_OUT <- resampled_shapes


all_outlines_OUT_centered <- Momocs::coo_centre(all_outlines_OUT) # center
all_outlines_OUT_centered_scaled <- Momocs::coo_scale(all_outlines_OUT_centered) # scale

all_outlines_OUT_centered_scaled <- Momocs::coo_slidedirection(all_outlines_OUT_centered_scaled, 
                                                               direction = "up")

# stack(all_outlines_OUT_centered_scaled)
# Momocs::panel(all_outlines_OUT_centered_scaled, fac = "Dataset")


all_outlines_OUT_centered_scaled_harmonics <- Momocs::calibrate_harmonicpower_efourier(all_outlines_OUT_centered_scaled,plot = F)
all_outlines_OUT_centered_scaled_efourier <- Momocs::efourier(all_outlines_OUT_centered_scaled,
                                                              nb.h = as.matrix(all_outlines_OUT_centered_scaled_harmonics[["minh"]])[[4,1]], norm = F) 
all_outlines_OUT_centered_scaled_PCA <- Momocs::PCA(all_outlines_OUT_centered_scaled_efourier) # PCA on Coe objects, using prcomp.

# plot(all_outlines_OUT_centered_scaled_PCA,
#      xax = 1,
#      yax = 2,
#      pos.shp = "XY",
#      rug = FALSE,
#      zoom = 1,
#      lwd.shp = 1,
#      size.shp = 0.5,
#      amp.shp = 0.5,
#      title = "",
#      chull = F,
#      "Dataset") # is NOT the 95% confidence ellipse, but just a convex hull around all group members
# 
# 
# minimum_no_of_pcs <- Momocs::scree_min(all_outlines_OUT_centered_scaled_PCA,
#                                        prop = 0.95) 
# minimum_no_of_pcs
# 
# Momocs::scree_plot(all_outlines_OUT_centered_scaled_PCA,
#                    nax = 1:(minimum_no_of_pcs+2))
# 
# minimum_no_of_pcs <- 4
# 
# # PCA shape variation
# gg <- Momocs::PCcontrib(all_outlines_OUT_centered_scaled_PCA,
#                         nax = 1:4, # sind alle minimum_no_of_pcs_petrik PCs "signifikant"/beschreiben alle PCs eine morphologische Eigenschaft? 
#                         sd.r = c(-2,-1,0,1,2)) 
# pc_contrib_plot <- gg$gg +
#   theme_bw() +
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank())
# 
# PCA_df_names <- as.data.frame(all_outlines_OUT_centered_scaled_PCA$x[,c(1:minimum_no_of_pcs)])
# PCA_df_names$ID <- row.names(PCA_df_names)
# 
# all_outlines_OUT_PCA_df <- dplyr::left_join(PCA_df_names, all_outlines_OUT$fac, by = "ID")
# 
# a <- ggplot(data = all_outlines_OUT_PCA_df, aes(fill = Dataset)) +
#   geom_point(pch = 21, size = 3) +
#   coord_fixed(ratio =1) +
#   theme_classic() +
#   theme(legend.position = "right",
#         axis.title = element_text(size = 14),
#         axis.text = element_text(size = 12),
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 16)) +
#   geom_hline(yintercept=0, linetype="dashed", alpha = 0.5) + 
#   geom_vline(xintercept=0, linetype="dashed", alpha = 0.5) +
#   ggthemes::scale_fill_colorblind()
# 
# 
# b <- a + aes(x = PC1, y = PC2) +
#   xlab(paste0("PC1 (", round(all_outlines_OUT_centered_scaled_PCA$eig[1]*100, digits = 0), "%)")) +
#   ylab(paste0("PC2 (", round(all_outlines_OUT_centered_scaled_PCA$eig[2]*100, digits = 0), "%)")) +
#   guides(color = FALSE, fill = FALSE)
# 
# c <- a + aes(x = PC1, y = PC3) +
#   xlab(paste0("PC1 (", round(all_outlines_OUT_centered_scaled_PCA$eig[1]*100, digits = 0), "%)")) +
#   ylab(paste0("PC3 (", round(all_outlines_OUT_centered_scaled_PCA$eig[3]*100, digits = 0), "%)")) + 
#   labs(fill = "Dataset")  +
#   guides(color = FALSE) 
# 
# d <- a + aes(x = PC1, y = PC4) +
#   xlab(paste0("PC1 (", round(all_outlines_OUT_centered_scaled_PCA$eig[1]*100, digits = 0), "%)")) +
#   ylab(paste0("PC4 (", round(all_outlines_OUT_centered_scaled_PCA$eig[4]*100, digits = 0), "%)")) +
#   guides(color = FALSE, fill = FALSE)
# 
# 
# # cowplot::plot_grid(b,c,d,pc_contrib_plot)


#################################
####### OUTLIERS ###########
#################################
outliers_db <- fpc::dbscan(all_outlines_OUT_centered_scaled_PCA$x, eps = 0.25, MinPts = 3)
plot(outliers_db, all_outlines_OUT_centered_scaled_PCA$x, main = "DBSCAN", frame = FALSE)


outliers_cluster <- data.frame(name = row.names(all_outlines_OUT_centered_scaled_PCA$x), 
                               value = outliers_db$cluster, 
                               row.names = NULL)
table(outliers_cluster$value)

outlier_names <- outliers_cluster[outliers_cluster$value == 0,]$name

outlier_shapes <- Momocs::filter(all_outlines_OUT_centered_scaled,
                                 ID %in% outlier_names)
Momocs::panel(outlier_shapes)

all_outlines_OUT_centered_scaled_no_outliers <- Momocs::slice(all_outlines_OUT_centered_scaled, 
                                                              -match(outlier_names, all_outlines_OUT_centered_scaled$fac$ID))

all_outlines_OUT_centered_scaled_harmonics <- Momocs::calibrate_harmonicpower_efourier(all_outlines_OUT_centered_scaled_no_outliers,plot = F)
all_outlines_OUT_centered_scaled_efourier <- Momocs::efourier(all_outlines_OUT_centered_scaled_no_outliers,
                                                              nb.h = as.matrix(all_outlines_OUT_centered_scaled_harmonics[["minh"]])[[4,1]], norm = F) 
all_outlines_OUT_centered_scaled_PCA <- Momocs::PCA(all_outlines_OUT_centered_scaled_efourier) # PCA on Coe objects, using prcomp.

minimum_no_of_pcs <- Momocs::scree_min(all_outlines_OUT_centered_scaled_PCA,
                                       prop = 0.95) 
minimum_no_of_pcs


screeplot_PCA <- Momocs::scree_plot(all_outlines_OUT_centered_scaled_PCA,
                   nax = 1:(minimum_no_of_pcs+1)) + 
  theme_bw() +
  xlab("Components") +
  theme(legend.position = "right",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))
screeplot_PCA

# PCA shape variation
gg <- Momocs::PCcontrib(all_outlines_OUT_centered_scaled_PCA,
                        nax = 1:4, # sind alle minimum_no_of_pcs_petrik PCs "signifikant"/beschreiben alle PCs eine morphologische Eigenschaft? 
                        sd.r = c(-2,-1,0,1,2)) 
pc_contrib_plot <- gg$gg +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

all_outlines_OUT_PCA_df <- as.data.frame(all_outlines_OUT_centered_scaled_PCA$x)
all_outlines_OUT_PCA_df$ID <- all_outlines_OUT_centered_scaled_PCA$fac$ID

all_outlines_OUT_PCA_df <- dplyr::left_join(all_outlines_OUT_PCA_df, all_outlines_OUT_centered_scaled_PCA$fac, by = "ID")

a <- ggplot(data = all_outlines_OUT_PCA_df, aes(fill = Dataset)) +
  geom_point(pch = 21, size = 3) +
  coord_fixed(ratio =1) +
  theme_classic() +
  theme(legend.position = "right",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)) +
  geom_hline(yintercept=0, linetype="dashed", alpha = 0.5) + 
  geom_vline(xintercept=0, linetype="dashed", alpha = 0.5) 


b <- a + aes(x = PC1, y = PC2) +
  xlab(paste0("PC1 (", round(all_outlines_OUT_centered_scaled_PCA$eig[1]*100, digits = 0), "%)")) +
  ylab(paste0("PC2 (", round(all_outlines_OUT_centered_scaled_PCA$eig[2]*100, digits = 0), "%)")) +
  guides(color = FALSE, fill = FALSE)

c <- a + aes(x = PC1, y = PC3) +
  xlab(paste0("PC1 (", round(all_outlines_OUT_centered_scaled_PCA$eig[1]*100, digits = 0), "%)")) +
  ylab(paste0("PC3 (", round(all_outlines_OUT_centered_scaled_PCA$eig[3]*100, digits = 0), "%)")) + 
  labs(fill = "Dataset")  +
  guides(color = FALSE, fill = FALSE) 

d <- a + aes(x = PC1, y = PC4) +
  xlab(paste0("PC1 (", round(all_outlines_OUT_centered_scaled_PCA$eig[1]*100, digits = 0), "%)")) +
  ylab(paste0("PC4 (", round(all_outlines_OUT_centered_scaled_PCA$eig[4]*100, digits = 0), "%)")) +
  guides(color = FALSE, fill = FALSE)


cowplot::plot_grid(b,c,d,pc_contrib_plot)

######################################################
# custom bins based on region
rownames_DATASETS <- list()
for(i in unique(all_outlines_OUT_PCA_df$Dataset)){
  rownames_DATASETS[[i]] <- as.character(subset(all_outlines_OUT_centered_scaled_PCA$fac, Dataset == i)$ID)
}

library(dispRity)
TS_subsets <- custom.subsets(all_outlines_OUT_centered_scaled_PCA$x, 
                             group = rownames_DATASETS)

TS_boot <- boot.matrix(TS_subsets, bootstraps = 1000)
TS_disp <- dispRity(TS_boot, metric = c(sum, variances))
summary(TS_disp)

summary(TS_disp)[,c(1,3)] # mean disparity 

# Wilcox.test
test.dispRity(TS_disp, 
              test = wilcox.test, 
              comparisons = "pairwise",
              correction = "bonferroni")
# PERMANOVA
test.dispRity(TS_disp, 
              test = adonis.dispRity, 
              comparisons = "pairwise",
              correction = "bonferroni")


TS_names <- names(TS_disp$disparity)
disparity_df_list <- list()
for(i in TS_names){

  if (i == "Palaeoindian" || i == "Bell Beaker CZ" || i == "Tanged points"){
    disparity_df_list[[i]] <- data.frame(Subset = paste0(i, 
                                                          "\n(n=",nrow(TS_disp$subsets[[i]]$elements),")"),
                                         disparity = as.vector(TS_disp$disparity[[i]][[2]]),
                                         nelements = nrow(TS_disp$subsets[[i]]$elements),
                                         TS = i)
  } else {
    disparity_df_list[[i]] <- data.frame(Subset = paste0(i, 
                                                          " (n=",nrow(TS_disp$subsets[[i]]$elements),")"),
                                         disparity = as.vector(TS_disp$disparity[[i]][[2]]),
                                         nelements = nrow(TS_disp$subsets[[i]]$elements),
                                         TS = i)
  }
  
}
disparity_df_TSdiscrete_armatureOutlines_perTShard <- do.call(rbind.data.frame, disparity_df_list)

disparity_TSdiscrete_armatureOutlines_ggplot_perTShard <- ggplot(data = disparity_df_TSdiscrete_armatureOutlines_perTShard, aes(x = Subset, y = disparity)) +
  geom_violin() + 
  geom_boxplot(notch = T, width = 0.1, aes(fill = TS)) +
  theme_bw() +
  ggtitle(NULL) +
  xlab("") + 
  ylab("Disparity (sum of variances)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text=element_text(size=14), #,face="bold"
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
        axis.title=element_text(size=14)) +
  # ggthemes::scale_fill_colorblind() +
  guides(color = FALSE, fill = FALSE)
disparity_TSdiscrete_armatureOutlines_ggplot_perTShard


whole_grid <- cowplot::plot_grid(b, screeplot_PCA, c, pc_contrib_plot, d, disparity_TSdiscrete_armatureOutlines_ggplot_perTShard,
                                 align = "hv", ncol = 2, labels = "AUTO", axis = "l")
whole_grid

height <- 380
width <- 285
# ggsave(whole_grid,
#        filename = file.path(".", "3_output", "whole_grid_resampled_by_region.svg"),
#        height  = height,
#        width = width,
#        units = "mm")
# ggsave(whole_grid,
#        filename = file.path(".", "3_output", "whole_grid_resampled_by_region.png"),
#        height  = height,
#        width = width,
#        units = "mm")




