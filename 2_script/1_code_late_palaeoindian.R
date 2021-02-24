library(readr)
library(ggplot2)
library(ggpubr)
library(Momocs)
library(caret)
library(vegan)
library(parallel)
library(MVN)
library(psych)



rm(list=ls())

output_folder <- file.path("./3_output/late_palaeoindian_bifac_points")
dir.create(output_folder,
           recursive = T)


# load outlines
outlines_combined_buchanan <- readRDS(file.path("./1_data/outlines_combined_goshen_plainview.RDS"))

# rename artefact names to match spreadsheet
buchanan_names <- list()
for (i in 1:length(outlines_combined_buchanan$coo)){
  buchanan_names[i] <- strsplit(names(outlines_combined_buchanan$coo), split = "rotated-")[[i]][2]
}
names(outlines_combined_buchanan$coo) <- as.vector(do.call(rbind, buchanan_names))

# load spreadsheet
buchanan_supp <- as.data.frame(readr::read_csv(file.path("./1_data/Buchanan_et_al_AA_Supp_Mats_edited.csv"),
                                               col_types = cols(prev_assigned_type = col_factor(levels = c("Goshen", "Plainview")))))

# add spreadsheet data to outlines
outlines_combined_buchanan <- Momocs::Out(outlines_combined_buchanan$coo,
                                          fac = buchanan_supp)


# GMM procedures
outlines_combined_buchanan_centered <- Momocs::coo_centre(outlines_combined_buchanan) # center
outlines_combined_buchanan_centered_scaled <- Momocs::coo_scale(outlines_combined_buchanan_centered) # scale
outlines_combined_buchanan_centered_scaled <- Momocs::coo_slidedirection(outlines_combined_buchanan_centered_scaled,
                                                                         direction = "up") 
outlines_combined_centered_scaled_smoothed <- outlines_combined_buchanan_centered_scaled


# # outline inspection
# stack(outlines_combined_buchanan_centered_scaled)

Momocs::panel(outlines_combined_buchanan_centered_scaled, fac = "prev_assigned_type")
# Goshen-type outlines
Momocs::panel(Momocs::slice(outlines_combined_buchanan_centered_scaled, prev_assigned_type == "Goshen"),
                main = paste("Goshen, n = ", nrow(filter(outlines_combined_buchanan_centered_scaled$fac, prev_assigned_type == "Goshen"))))
# Plainview-type outlines
Momocs::panel(Momocs::slice(outlines_combined_buchanan_centered_scaled, prev_assigned_type == "Plainview"),
                main = paste("Plainview, n = ", nrow(filter(outlines_combined_buchanan_centered_scaled$fac, prev_assigned_type == "Plainview"))))

  

# harmonic calibration. Estimates the number of harmonics required for the Fourier methods implemented in Momocs. This is the only step in this section that produces data we need in the subsequent step.
outlines_combined_buchanan_centered_scaled_harmonics <- Momocs::calibrate_harmonicpower_efourier(outlines_combined_buchanan_centered_scaled, 
                                                                                                 plot = F)
# outlines_combined_buchanan_centered_scaled_harmonics
# efourier
outlines_combined_buchanan_centered_scaled_efourier <- Momocs::efourier(outlines_combined_buchanan_centered_scaled,
                                                                        nb.h = as.matrix(outlines_combined_buchanan_centered_scaled_harmonics[["minh"]])[[4,1]], # choses number of harmonics for 99.9%
                                                                        norm = F) 
# PCA
outlines_combined_buchanan_centered_scaled_PCA <- Momocs::PCA(outlines_combined_buchanan_centered_scaled_efourier) # PCA on Coe objects, using prcomp.

minimum_no_of_pcs_buchanan <- Momocs::scree_min(outlines_combined_buchanan_centered_scaled_PCA,
                                                prop = 0.95) 
minimum_no_of_pcs_buchanan # to describe 95% of the data's variability

buchanan_screeplot <- Momocs::scree_plot(outlines_combined_buchanan_centered_scaled_PCA,
                   nax = 1:8) + theme_bw()
ggsave(buchanan_screeplot,
       filename = file.path("./3_output/late_palaeoindian_bifac_points/buchanan_screeplot.svg"),
       width = 7,
       height = 5,
       dpi = 320,
       units = "in",
       device = "svg")
ggsave(buchanan_screeplot,
       filename = file.path("./3_output/late_palaeoindian_bifac_points/buchanan_screeplot.png"),
       width = 7,
       height = 5,
       dpi = 320,
       units = "in",
       device = "png")


# PCA shape variation
buchanan_outline_analysis_shp_variation <- Momocs::PCcontrib(outlines_combined_buchanan_centered_scaled_PCA,
                                                             nax = 1:minimum_no_of_pcs_buchanan, # sind alle minimum_no_of_pcs_buchanan PCs "signifikant"/beschreiben alle PCs eine morphologische Eigenschaft? 
                                                             sd.r = c(-2,-1,0,1,2)) 


buchanan_outline_analysis_shp_variation_plot <- buchanan_outline_analysis_shp_variation$gg + 
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave(buchanan_outline_analysis_shp_variation_plot,
       filename = file.path("./3_output/late_palaeoindian_bifac_points/buchanan_outline_analysis_shp_variation_plot.svg"),
         width = 7,
       height = 5,
       dpi = 320,
       units = "in",
       device = "svg")
ggsave(buchanan_outline_analysis_shp_variation_plot,
       filename = file.path("./3_output/late_palaeoindian_bifac_points/buchanan_outline_analysis_shp_variation_plot.png"),
       width = 7,
       height = 5,
       dpi = 320,
       units = "in",
       device = "png")



# PCA plot
buchanan_PCA_data <- data.frame(cbind(outlines_combined_buchanan_centered_scaled_PCA$x[,1:minimum_no_of_pcs_buchanan], 
                                      buchanan_supp$prev_assigned_type))
names(buchanan_PCA_data)[ncol(buchanan_PCA_data)] <- "prev_assigned_type"

buchanan_PCA_data_names <- buchanan_PCA_data
buchanan_PCA_data_names$prev_assigned_type[buchanan_PCA_data_names$prev_assigned_type == "1"] <- "Goshen"
buchanan_PCA_data_names$prev_assigned_type[buchanan_PCA_data_names$prev_assigned_type == "2"] <- "Plainview"

PCA_plot_GoshenPlainview <- ggplot(data = buchanan_PCA_data_names, aes(x = PC1, y = PC2, fill = prev_assigned_type)) +
  geom_point(pch = 21, size = 3) +
  stat_ellipse(aes(color = prev_assigned_type)) + # https://stats.stackexchange.com/a/38506
  coord_fixed(ratio =1) +
  theme_classic() +
  xlab(paste0("PC1 (", round(outlines_combined_buchanan_centered_scaled_PCA$eig[1]*100, digits = 0), "%)")) +
  ylab(paste0("PC2 (", round(outlines_combined_buchanan_centered_scaled_PCA$eig[2]*100, digits = 0), "%)")) +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)) +
  labs(fill = "Type") +
  guides(color = FALSE) + 
  geom_hline(yintercept=0, linetype="dashed", alpha = 0.5) + 
  geom_vline(xintercept=0, linetype="dashed", alpha = 0.5)
PCA_plot_GoshenPlainview

ggsave(PCA_plot_GoshenPlainview,
       filename = file.path("./3_output/late_palaeoindian_bifac_points/PCA_plot_GoshenPlainview.svg"),
         width = 7,
       height = 5,
       dpi = 320,
       units = "in",
       device = "svg")
ggsave(PCA_plot_GoshenPlainview,
       filename = file.path("./3_output/late_palaeoindian_bifac_points/PCA_plot_GoshenPlainview.png"),
       width = 7,
       height = 5,
       dpi = 320,
       units = "in",
       device = "png")


# Discriminant function analysis is broken into a 2-step process: (1) testing significance of a set of discriminant functions, and; (2) classification. The first step is computationally identical to MANOVA. There is a matrix of total variances and covariances; likewise, there is a matrix of pooled within-group variances and covariances. The two matrices are compared via multivariate F tests in order to determine whether or not there are any significant differences (with regard to all variables) between groups. One first performs the multivariate test, __and, if statistically significant__, proceeds to see which of the variables have significantly different means across the groups. Once group means are found to be statistically significant, classification of variables is undertaken [@poulsen_discriminant_nodate].
# 0. does the data follow a normal distribution?
buchanan_PCA_data_GOSHEN <- subset(buchanan_PCA_data, prev_assigned_type == 1)
buchanan_PCA_data_PLAINVIEW <- subset(buchanan_PCA_data, prev_assigned_type == 2)

psych::multi.hist(buchanan_PCA_data_GOSHEN[,c(1:minimum_no_of_pcs_buchanan)])
MVN::mvn(data = as.matrix(buchanan_PCA_data_GOSHEN[,c(1:minimum_no_of_pcs_buchanan)]), 
         multivariatePlot = "qq")

psych::multi.hist(buchanan_PCA_data_PLAINVIEW[,c(1:minimum_no_of_pcs_buchanan)])
MVN::mvn(data = as.matrix(buchanan_PCA_data_PLAINVIEW[,c(1:minimum_no_of_pcs_buchanan)]), 
         multivariatePlot = "qq")

# 1. MANOVA with Lawley-Hotelling statistic, also known as Hotelling’s generalized T^2-statistic
# $H_{0}$: The data are drawn from populations with the same multivariate means.
# However, the test assumes a multivariate normal distribution, which is violated in this case (see above)
# MANOVA is not robust against unequal sample sizes! If it it the smaller sample that produces larger variances and covariances, the probability values will be liberal and so significant differences should be treated with caution (although non-significant effects can be trusted) [@field_discovering_2012, 718]
regular_manova <- summary(manova(as.matrix(buchanan_PCA_data[,c(1:minimum_no_of_pcs_buchanan)])~as.factor(buchanan_supp$prev_assigned_type)),
                          test="Hotelling")
regular_manova
if(round(regular_manova$stats[,"Pr(>F)"][[1]], digits = 3)<0.05){result_manova <- "had to be rejected"} else {result_manova <- "could not be rejected"}
print(paste0("H0 ", result_manova))
# 1.b) because of non-normal multivariate distribution, we conduct a permutation test of (squared) Mahalanobis distance, as described in @hammer_paleontological_2006, 66f.
# Distance calculation: Using the Mahalanobis distance between the groups centroids from the principal components.
# # Mahalanobis distance between the groups centroids from the principal components
buchanan_mahal_dist_GOSHEN_PLAINVIEW <- asbio::D.sq(buchanan_PCA_data_GOSHEN[,c(1:minimum_no_of_pcs_buchanan)], buchanan_PCA_data_PLAINVIEW[,c(1:minimum_no_of_pcs_buchanan)])
buchanan_mahal_dist_GOSHEN_PLAINVIEW$D.sq # Mahalanobis distance between the groups centroids based on the Principal Components data
round(buchanan_mahal_dist_GOSHEN_PLAINVIEW$D.sq, digits = 3)

# ## LDA based on PC-scores
library(splitstackshape)
set.seed(1)
stratified_sample <- splitstackshape::stratified(buchanan_PCA_data_names, 
                                                 group = "prev_assigned_type", 
                                                 size = 0.75, 
                                                 replace = F,
                                                 bothSets = T)

training <- as.data.frame(stratified_sample$SAMP1)
training$prev_assigned_type <- as.factor(training$prev_assigned_type)
training_wo_cluster <- training
training_wo_cluster$prev_assigned_type <- NULL

testing  <- as.data.frame(stratified_sample$SAMP2)
testing$prev_assigned_type <- as.factor(testing$prev_assigned_type)

# cross-validated misclassification matrix based on PCA data
library(caret)
buchanan_lda_CV <- caret::train(prev_assigned_type ~.,
                                data=training, 
                                method="lda",
                                trControl = trainControl(method = "cv"))


predicted_testing <- predict(buchanan_lda_CV, testing[,c(1:(ncol(testing)-1))])
summary_buchanan_lda_CV <- caret::confusionMatrix(predicted_testing, testing$prev_assigned_type,
                                                  mode = "prec_recall",
                                                  positive = c("Goshen"))
confusion_matrix_buchanan_lda_CV_df <- as.data.frame.matrix(summary_buchanan_lda_CV$table)
confusion_matrix_buchanan_lda_CV_df

buchanan_lda_variate_contribution <- data.frame(PC = row.names(as.data.frame(buchanan_lda_CV$finalModel$scaling)),
                                                round(as.data.frame(buchanan_lda_CV$finalModel$scaling), digits = 2),
                                                LD1_abs = abs(round(buchanan_lda_CV$finalModel$scaling[,1], digits = 2)))


buchanan_lda_variate_contribution_reordered <- buchanan_lda_variate_contribution[order(buchanan_lda_variate_contribution$LD1_abs, decreasing = T),]
buchanan_lda_variate_contribution_reordered


# PERMANOVA on raw PC-scores; uses a permutation test with pseudo-F-ratios
n_permut <- 10000
permut_manova_buchanan_pca_euclid <- vegan::adonis(as.matrix(buchanan_PCA_data[,c(1:minimum_no_of_pcs_buchanan)]) ~ buchanan_PCA_data$prev_assigned_type,
                                                   method = "euclidean", # method to calculate pairwise distances
                                                   permutations = n_permut,
                                                   parallel = parallel::detectCores()-1)
permut_manova_buchanan_pca_euclid
if (round(permut_manova_buchanan_pca_euclid$aov.tab$"Pr(>F)"[1], digits = 3) < 0.05){result_permut_manova_buchanan_pca_euclid <- "had to be rejected, meaning the means of all groups are not equal"} else {result_permut_manova_buchanan_pca_euclid <- "could not be rejected, meaning the means of all groups are equal"}
result_permut_manova_buchanan_pca_euclid
####




# # measure of relative shape distances as proposed by Klingenberg and Monteiro (2005) [@klingenberg_distances_2005]
# # the distance measure provides a "scalar measure of the relative extent of shape differences, while taking into account that variation may not be isotropic. This distance measure is based on the idea of one-sample standard distance (Flury and Redwyl, 1986; FLury, 1997), which is equivalent to the one-sample version of the Mahalanobis distance (Mardia et al., 1979: 31)."

# #1. carry out a principal components analysis
# # done: 
# outlines_combined_buchanan_centered_scaled_PCA$x[,1:minimum_no_of_pcs_buchanan]

# 2. standardize the scores: subtract the mean and divide by the standard deviation
standardization_fun <- function(obs = X, mean = mean, sd = sd){
  standarized_obs <- (obs-mean)/sd
  return(standarized_obs)
}

pc_scores_mean_list <- list()
pc_scores_sd_list <- list()
pc_scores_standardized_and_squared_list <- list()

for (col_index in 1:minimum_no_of_pcs_buchanan){
  pc_scores_mean_list[[col_index]] <- mean(outlines_combined_buchanan_centered_scaled_PCA$x[,col_index])
  mean <- pc_scores_mean_list[[col_index]]
  pc_scores_sd_list[[col_index]] <- sd(outlines_combined_buchanan_centered_scaled_PCA$x[,col_index])
  sd <- pc_scores_sd_list[[col_index]]
  
  pc_scores_standardized_and_squared_list[[col_index]] <- (standardization_fun(as.data.frame(outlines_combined_buchanan_centered_scaled_PCA$x[,col_index]),
                                                                               mean = mean,
                                                                               sd = sd))^2
}
standardized_and_squared_scores_PCA_buchanan <- do.call(cbind.data.frame, pc_scores_standardized_and_squared_list)

for (col_index in 1:minimum_no_of_pcs_buchanan){
  names(standardized_and_squared_scores_PCA_buchanan)[col_index] <- paste0("PC", col_index)
}

# 3. sum up the squares of these standardized principal component scores for every observation
standardized_and_squared_scores_summed_PCA_buchanan <- as.data.frame(apply(standardized_and_squared_scores_PCA_buchanan, 1, sum))

# 4. compute the square root of the resulting sum for each observation
sqrt_of_standardized_and_squared_scores_summed_PCA_buchanan <- as.data.frame(apply(standardized_and_squared_scores_summed_PCA_buchanan, 1, sqrt))
names(sqrt_of_standardized_and_squared_scores_summed_PCA_buchanan)[1] <- "klingenberg_shp_dist"
sqrt_of_standardized_and_squared_scores_summed_PCA_buchanan$prev_assigned_type <- buchanan_supp$prev_assigned_type


# 4.5 testing whether both samples are from the same distribution
# A) permutated non-paramteric t-test (Mann-Whitney test) from https://rpubs.com/deleeuw/9319
t.np <- function(x, y, p = 0.95, nperm = nperm) {
  n1 <- length(x)
  n2 <- length(y)
  t0 <- abs(t.stat(x, y))
  nn <- 1:(n1 + n2)
  z <- c(x, y)
  tt <- NULL
  for (i in 1:nperm) {
    ind <- sample(nn, n1)
    xp <- z[ind]
    yp <- z[-ind]
    tt <- c(tt, abs(t.stat(xp, yp)))
  }
  pp <- edf(t0, tt)
  return(list(sig = pp > p, pp = pp))
}
# function to compute the t-statistic for two-sided testing
t.stat <- function(x, y) {
  n1 <- length(x)
  n2 <- length(y)
  mx <- mean(x)
  my <- mean(y)
  nn <- n1 + n2 - 2
  nm <- (1/n1) + (1/n2)
  sp <- sqrt(((n1 - 1) * var(x) + (n2 - 1) * var(y)) * nm/nn)
  tt <- (mx - my)/sp
  return(tt)
}
# For our non-parametric computations we need a simple auxilary to compute the empirical distribution function.
edf <- function(z, x) {
  sapply(z, function(y) length(which(x < y))/length(x))
}

# permutated non-paramteric t-test (Mann-Whitney test) from https://rpubs.com/deleeuw/9319
# t.np(subset(sqrt_of_standardized_and_squared_scores_summed_PCA_buchanan, prev_assigned_type == "Goshen")$klingenberg_shp_dist,
#      subset(sqrt_of_standardized_and_squared_scores_summed_PCA_buchanan, prev_assigned_type == "Plainview")$klingenberg_shp_dist,
#      nperm = 10000)

# B) Mann-Whitney-U test as non-paramtetric test
buchanan_relat_shp_diff_mannwhitney <- wilcox.test(sqrt_of_standardized_and_squared_scores_summed_PCA_buchanan$klingenberg_shp_dist~sqrt_of_standardized_and_squared_scores_summed_PCA_buchanan$prev_assigned_type)
# buchanan_relat_shp_diff_mannwhitney

# C) permutated independece test
# coin::independence_test(sqrt_of_standardized_and_squared_scores_summed_PCA_buchanan$klingenberg_shp_dist~sqrt_of_standardized_and_squared_scores_summed_PCA_buchanan$prev_assigned_type,
#                   data = sqrt_of_standardized_and_squared_scores_summed_PCA_buchanan,
#                   distribution = approximate(nresample = 2000))


# 5. histogram-plot 
buchanan_klingenberg_hist <- ggplot(sqrt_of_standardized_and_squared_scores_summed_PCA_buchanan) +
  geom_histogram(aes(x=klingenberg_shp_dist, fill=prev_assigned_type), 
                 colour="grey50",
                 position="dodge") +
  labs(x="Relative shape distance", 
       y="Count") +
  labs(fill="Type") +
  annotate("label", 
           x = 4, 
           y = 12, 
           label = paste0("Mann-Whitney-U test:\nW = ", buchanan_relat_shp_diff_mannwhitney$statistic[[1]],", p = ", round(buchanan_relat_shp_diff_mannwhitney$p.value, digits = 3))) + 
  theme_bw()

ggsave(filename = file.path("./3_output/late_palaeoindian_bifac_points/buchanan_pca_rel_shp_klingenberg_hist.svg"),
       plot = buchanan_klingenberg_hist,
       device = "svg",
       width = 10,
       height = 10,
       dpi = "retina")
ggsave(filename = file.path("./3_output/late_palaeoindian_bifac_points/buchanan_pca_rel_shp_klingenberg_hist.png"),
       plot = buchanan_klingenberg_hist,
       device = "png",
       width = 10,
       height = 10,
       dpi = "retina")


buchanan_klingenberg_hist_descr <- buchanan_klingenberg_hist +
  labs(caption = paste0("Distribution of measure of relative shape distance (Klingenberg and Monteiro, 2005) between Goshen (n=", nrow(buchanan_PCA_data_GOSHEN), ") and Plainview (n=", nrow(buchanan_PCA_data_PLAINVIEW), ") type outlines."))

ggsave(filename = file.path("./3_output/late_palaeoindian_bifac_points/buchanan_pca_rel_shp_klingenberg_hist_descr.svg"),
       plot = buchanan_klingenberg_hist_descr,
       device = "svg",
       width = 10,
       height = 10,
       dpi = "retina")
ggsave(filename = file.path("./3_output/late_palaeoindian_bifac_points/buchanan_pca_rel_shp_klingenberg_hist_descr.png"),
       plot = buchanan_klingenberg_hist_descr,
       device = "png",
       width = 10,
       height = 10,
       dpi = "retina")
