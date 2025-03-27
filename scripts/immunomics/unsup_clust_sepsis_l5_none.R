rm(list = ls())
#remove.packages("plotUtils")
#detach("package:plotUtils")
if (!require("plotUtils", quietly = T)){
        devtools::install_github("guisantagui/plotUtils", upgrade = "never")
}
plotUtils::plotPCA
library(plotUtils)
library(ggplot2)
subgroup_data_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/data/TIPS patient classification_sb_cs27022025.xlsx"
group_data_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/data/TIPS-Data_analysis/TIPS-Results_collection/TIPS_data_collection_2025-01-11/grouping_master_data_2025-01-10.RDS"
immu_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/data/Immunomics_master_data_2025-01-10.RData"
outDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/results/immunomics/"
outDir <- add_slash_if_not(outDir)

create_dir_if_not(outDir)

load(immu_f)
group_data <- readRDS(group_data_f)
subgroup_data <- as.data.frame(readxl::read_xlsx(subgroup_data_f))

comb_df <- cbind.data.frame(group_data,
                            subgroup_data[match(group_data$patient_ID,
                                                subgroup_data$`Study id`),
                                          c("subgroup", "regrouping")])
comb_df <- comb_df[!is.na(comb_df$group), ]
write.csv(comb_df, file = sprintf("%sgrouping_subgroups.csv", outDir))
#View(immunomics_counts_level2)



not_in_groupDat <- unique(immunomics_counts_level1$patient_ID[!immunomics_counts_level1$patient_ID %in% group_data$patient_ID])

parse_immun <- function(immun, patient_info = NULL, type = "counts"){
        #immun <- immunomics_counts_level1
        #patient_info <- group_data
        immun_comm <- data.frame(immun[, 1:grep("dalama", colnames(immun))])
        if (type == "counts"){
                immun_dat <- data.frame(immun[, (grep("dalama",
                                                      colnames(immun)) + 2):ncol(immun)])
                immun_dat_norm <- (immun_dat + 1)/(immun$allcounts + 1)
        }else if (type == "states"){
                immun_dat <- data.frame(immun[, (grep("dalama",
                                                      colnames(immun)) + 1):ncol(immun)])
                immun_dat <- as.data.frame(apply(immun_dat, 2, as.numeric))
        }
        if (!is.null(patient_info)){
                immun_comm$group <- patient_info$group[match(immun_comm$patient_ID,
                                                             patient_info$patient_ID)]
        }
        rownames(immun_dat) <- immun_comm$immunomics_ID
        colnames(immun_comm) <- gsub("immunomics_ID",
                                     "sample",
                                     colnames(immun_comm))
        
        if (type == "counts"){
                rownames(immun_dat_norm) <- rownames(immun_dat)
                out <- list(sample_info = immun_comm,
                            data = immun_dat,
                            norm = immun_dat_norm)
        }else if (type == "states"){
                out <- list(sample_info = immun_comm,
                            data = immun_dat)
        }
        return(out)
}

do_log2 <- function(obj, slot = "data", add_base = 0){
        obj$log2 <- log2(obj[[slot]] + add_base)
        return(obj)
}

do_stand <- function(obj, slot = "data"){
        obj$stand <- stand(obj[[slot]])
        return(obj)
}

filt_samps <- function(obj, keep_var = NULL,
                       keep_val = NULL,
                       keep_samps = NULL){
        if (is.null(keep_var) & is.null(keep_val)){
                if (is.null(keep_samps)){
                        stop("If keep_var and keep_val are NULL a vector keep_samp with samples to keep needs to be given.",
                             call. = F)
                }
        }else if (xor(is.null(keep_var), is.null(keep_val))){
                stop("Both keep_var and keep_val need to be provided for them to be used",
                     call. = F)
        }else if (!is.null(keep_var) & !is.null(keep_val)){
                if (!is.null(keep_samps)){
                        warning("As keep_var and keep_val are provided, keep_samps won't be used.",
                                call. = F)
                }
                keep_samps <- obj$sample_info$sample[obj$sample_info[, keep_var] %in% keep_val]
        }
        obj$sample_info <- obj$sample_info[obj$sample_info$sample %in% keep_samps, ]
        slots <- names(obj)
        slots <- slots[slots != "sample_info"]
        for (s in slots){
                obj[[s]] <- obj[[s]][rownames(obj[[s]]) %in% obj$sample_info$sample, ]
        }
        return(obj)
}

# Creates a list of multi-scoreplots from a list of PCAs
multPlot_pca_list <- function(pca_lst,
                              samp_info_list,
                              col = NULL,
                              nComps = 5,
                              biplot = T,
                              topNFeats = 5,
                              point_size = 3){
        plot_list <- mapply(function(pca, si) plotUtils::doPCAMultiPlot(pca,
                                                                        samp_info = si,
                                                                        col = col,
                                                                        nComps = 5,
                                                                        biplot = biplot,
                                                                        topNFeats = topNFeats,
                                                                        point_size = point_size),
                            pca_lst,
                            samp_info_list,
                            SIMPLIFY = F)
        return(plot_list)
}

# Obtains a dataframe of the cumulative distribution functions obtained
# from the consensus matrixes of a consensus clustering object (generated
# with consensusClusterPlus R package)
getCDFValues <- function(consClustObjkt, resolution = .0001){
        cdfDF <- data.frame(matrix(nrow = 0,
                                   ncol = 3,
                                   dimnames = list(NULL,
                                                   c("k",
                                                     "consensus_index",
                                                     "CDF"))))
        areas <- c()
        for(i in 2:length(consClustObjkt)){
                k <- i
                vals <- c()
                xAxis <- seq(0, 1, resolution)
                for(j in xAxis){
                        consMat <- consClustObjkt[[k]]$consensusMatrix
                        upperTri <- consMat[upper.tri(consMat)]
                        val <- sum(upperTri <= j)/length(upperTri)
                        vals <- c(vals, val)
                }
                toBind <- data.frame(k = rep(k, length(vals)),
                                     consensus_index = xAxis,
                                     CDF = vals)
                cdfDF <- rbind.data.frame(cdfDF, toBind)
                # Calculate area under the curve 
                multVec <- c(1,
                             rep(2,
                                 nrow(toBind) - 2),
                             1)
                trapArea <- resolution/2*sum(toBind$CDF * multVec)
                areas <- c(areas, trapArea)
        }
        names(areas) <- as.character(unique(cdfDF$k))
        cdfDF$k <- factor(cdfDF$k)
        names(areas) <- levels(cdfDF$k)
        out <- list(cdfDF = cdfDF,
                    areas = areas)
        return(out)
}

# Gives the PACs (proportion of ambiguous clustering) given a DF of CDF
# curves.
getPAC <- function(cdfDF){
        pac_df <- data.frame(matrix(nrow = 0,
                                    ncol = 2,
                                    dimnames = list(NULL,
                                                    c("k", "PAC"))))
        for(k in levels(cdfDF$k)){
                cdf_k <- cdfDF[cdfDF$k == k, ]
                cdfk2 <- cdf_k$CDF[cdf_k$consensus_index == 0.9]
                cdfk1 <- cdf_k$CDF[cdf_k$consensus_index == 0.1]
                pac <- cdfk2 - cdfk1
                toBind <- data.frame(matrix(c(as.numeric(k), pac),
                                            nrow = 1,
                                            ncol = 2,
                                            dimnames = list(NULL, colnames(pac_df))))
                pac_df <- rbind.data.frame(pac_df, toBind)
        }
        return(pac_df)
}

# PLots the CDF curves, the deltas and the PACs given a consensus cluster
# object (generated with consensusClusterPlus R package)
getCDFCurves <- function(consClustObjkt, resolution = .0001){
        cdfOut <- getCDFValues(consClustObjkt = consClustObjkt,
                               resolution = resolution)
        cdfDF <- cdfOut$cdfDF
        areas <- cdfOut$areas
        
        CDF_plot <- ggplot(data = cdfDF,
                           mapping = aes(x = consensus_index,
                                         y = CDF,
                                         col = k)) +
                geom_line() +
                ylim(0, 1) +
                theme(axis.text.y = element_text(size=15),
                      axis.text.x = element_text(size=15),
                      panel.background = element_blank(),
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black"),
                      axis.line.y = element_line(colour = "black"),
                      panel.border = element_rect(colour = "black",
                                                  fill=NA,
                                                  size=1))
        # Obtain Relative area under CDF curve change plot
        relAreaVec <- c()
        for(i in seq_along(areas)){
                if(i == 1){
                        relArea <- areas[i]
                }else{
                        relArea <- (areas[i] - areas[i - 1])/areas[i - 1]
                }
                relAreaVec <- c(relAreaVec, relArea)
        }
        relAreaDF <- data.frame(k = as.numeric(names(areas)),
                                rel_area_under_CDF_change = relAreaVec)
        
        relAreaPlot <- ggplot(data = relAreaDF, mapping = aes(x = k,
                                                              y = rel_area_under_CDF_change)) +
                geom_line() +
                geom_point() +
                theme(axis.text.y = element_text(size=15),
                      axis.text.x = element_text(size=15),
                      panel.background = element_blank(),
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black"),
                      axis.line.y = element_line(colour = "black"),
                      panel.border = element_rect(colour = "black",
                                                  fill=NA,
                                                  size=1))
        # Obtain PAC plot
        pac_df <- getPAC(cdfDF)
        pacPlot <- ggplot(data = pac_df, mapping = aes(x = k,
                                                       y = PAC)) +
                geom_line() +
                geom_point() +
                theme(axis.text.y = element_text(size=15),
                      axis.text.x = element_text(size=15),
                      panel.background = element_blank(),
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black"),
                      axis.line.y = element_line(colour = "black"),
                      panel.border = element_rect(colour = "black",
                                                  fill=NA,
                                                  size=1))
        out <- list(CDF = CDF_plot,
                    rel_CDF_area_change = relAreaPlot,
                    PAC = pacPlot)
        return(out)
}

# Obtain a plot including the consensus tree and the consensus matrix for
# a given consensus clust object and a given K.
getConsMatDendPlot <- function(consClustObjkt, k){
        consTree <- consClustObjkt[[k]]$consensusTree
        
        consMatx <- consClustObjkt[[k]]$consensusMatrix
        dimnames(consMatx) <- list(names(consClustObjkt[[k]]$consensusClass),
                                   names(consClustObjkt[[k]]$consensusClass))
        
        
        rowOrder <- consTree$order[length(consTree$order):1]
        consMatx <- consMatx[rowOrder, ]
        
        consClass <- consClustObjkt[[k]]$consensusClass
        
        classCols <- rainbow(max(consClass))[consClass]
        
        
        pdf(file = paste0(plotUnsupDir, sprintf("%s_consMatDend_k%s.pdf", outName, k)),
            width = 6, height = 6)
        heatmap.2(consMatx, dendrogram = "column",
                  Colv = as.dendrogram(consTree),
                  Rowv = F,
                  trace = "none",
                  notecol = NULL,
                  labRow = FALSE,
                  key = F,
                  col = colorpanel(75, low = "white", "blue"),
                  breaks = 76,
                  ColSideColors = classCols,
                  margins = c(6.5, 4))
        
        legend("left", title = "Groups",
               legend = sort(unique(consClass)),
               fill = rainbow(max(consClass)))
        dev.off()
}


# pca_red_var. If set to a value between 0 and 1, the dimensionality of the
# input data will be reduced with PCA to keep the indicated amount of variance
# before fitting the isolation forest.
detectOLs_isfrst <- function(dat,
                             quant_thrs = 0.95,
                             group_wise = F,
                             samp_info = NULL,
                             groups_in = NULL,
                             pca_red_var = NULL){
        if (!is.null(pca_red_var)){
                if (class(pca_red_var) != "numeric" | length(pca_red_var) != 1){
                        stop("pca_red_var has to be a numeric value.")
                }
                if (!(pca_red_var <= 1 & pca_red_var >= 0)){
                        stop("pca_red_var has to be between 0 and 1.")
                }
                dat_pca <- prcomp(dat, scale. = F, center = F)
                pca_imp <- summary(dat_pca)$importance
                max_pc_idx <- which(pca_imp["Cumulative Proportion", ] >= pca_red_var)[1]
                dat <- dat_pca$x[, 1:max_pc_idx]
        }
        if (group_wise){
                if (is.null(samp_info) | is.null(groups_in)){
                        stop("If group_wise == T, samp_info DF and variable where the groupings are (groups_in) needs to be provided.",
                             call. = F)
                }
                if (!"sample" %in% colnames(samp_info)){
                        stop("samp_info doesn't have a 'sample' column",
                             call. = F)
                }
                if (!groups_in %in% colnames(samp_info)){
                        stop(sprintf("%s isn't a column of samp_info",
                                     groups_in),
                             call. = F)
                }
                groups <- unique(samp_info[, groups_in])
                ols <- c()
                as_df <- data.frame(matrix(ncol = 2, nrow = 0,
                                           dimnames = list(NULL, c("sample",
                                                                   "anomaly_score"))))
                for (g in groups){
                        i_f <- isolation.forest(dat[samp_info$sample[samp_info[, groups_in] == g], ],
                                                ndim = 1,
                                                sample_size = 36,
                                                ntrees = 100,
                                                prob_pick_pooled_gain = 1)
                        if_pred <- predict(i_f, dat[samp_info$sample[samp_info[, groups_in] == g], ])
                        to_bind <- data.frame(sample = names(if_pred),
                                              anomaly_score = if_pred)
                        thrshld <- quantile(if_pred, probs = quant_thrs)
                        ols_group <- names(if_pred)[if_pred >= thrshld]
                        ols <- c(ols, ols_group)
                        as_df = rbind.data.frame(as_df, to_bind)
                }
        }else{
                i_f <- isolation.forest(dat,
                                        ndim = 1,
                                        sample_size = 36,
                                        ntrees = 100,
                                        prob_pick_pooled_gain = 1)
                if_pred <- predict(i_f, dat)
                thrshld <- quantile(if_pred, probs = quant_thrs)
                ols <- names(if_pred)[if_pred >= thrshld]
                as_df <- data.frame(sample = names(if_pred),
                                    anomaly_score = if_pred)
        }
        as_df <- as_df[order(as_df$anomaly_score, decreasing = T), ]
        out <- list(outliers = ols,
                    anomaly_scores = as_df)
        return(out)
}

################################################################################
################################################################################
# Counts analysis                                                              #
################################################################################
################################################################################

imm_counts_list <- list(l1 = immunomics_counts_level1,
                        l2 = immunomics_counts_level2,
                        l3 = immunomics_counts_level3,
                        l4 = immunomics_counts_level4,
                        l5 = immunomics_counts_level5)

# Parse the files
imm_counts_list_parsed <- lapply(imm_counts_list,
                                 parse_immun,
                                 patient_info = group_data,
                                 type = "counts")

# Filter to keep samples we have
imm_counts_list_parsed <- lapply(imm_counts_list_parsed,
                                 function(x) filt_samps(x,
                                                        keep_samps = x$sample_info$sample[!is.na(x$sample_info$group)]))

imm_counts_list_parsed <- lapply(imm_counts_list_parsed,
                                 do_log2,
                                 slot = "norm")
imm_counts_list_parsed <- lapply(imm_counts_list_parsed,
                                 do_stand,
                                 slot = "log2")

plot(density(apply(cbind.data.frame(lapply(imm_counts_list_parsed,
                                           function(x) x$data)),
                   2,
                   var)))

plot(density(apply(cbind.data.frame(lapply(imm_counts_list_parsed,
                                           function(x) x$stand)),
                   2,
                   var)))

sum(apply(cbind.data.frame(lapply(imm_counts_list_parsed,
                                  function(x) x$data)), 2, var) == 0)
quantile(apply(cbind.data.frame(lapply(imm_counts_list_parsed,
                                       function(x) x$data)), 2, var),
         probs = seq(0, 1, 0.1))


imm_counts_list_parsed_none <- lapply(imm_counts_list_parsed,
                                      function(x) filt_samps(x,
                                                             keep_var = "treatment",
                                                             keep_val = "none"))
imm_counts_list_parsed_LPS <- lapply(imm_counts_list_parsed,
                                     function(x) filt_samps(x,
                                                            keep_var = "treatment",
                                                            keep_val = "LPS"))
imm_counts_list_parsed_R848 <- lapply(imm_counts_list_parsed,
                                      function(x) filt_samps(x,
                                                             keep_var = "treatment",
                                                             keep_val = "R848"))


# Stand again
imm_counts_list_parsed_none <- lapply(imm_counts_list_parsed_none,
                                      do_stand,
                                      slot = "log2")
imm_counts_list_parsed_LPS <- lapply(imm_counts_list_parsed_LPS,
                                     do_stand,
                                     slot = "log2")
imm_counts_list_parsed_R848 <- lapply(imm_counts_list_parsed_R848,
                                      do_stand,
                                      slot = "log2")


# Perform consensus clustering on sepsis groups to infer subtypes
library(ConsensusClusterPlus)
library(M3C)

# Use Level 5 none, as its the most informative set.

l5_none_sepsis <- filt_samps(imm_counts_list_parsed_none$l5,
                             keep_var = "group",
                             keep_val = "sepsis")

l5_none_sepsis <- do_stand(l5_none_sepsis, slot = "log2")

consClust <- ConsensusClusterPlus(as.matrix(t(l5_none_sepsis$stand)),
                                  pItem = 0.8,
                                  pFeature = 1,
                                  clusterAlg = "hc",
                                  distance = "spearman",
                                  innerLinkage = "ward.D2",
                                  finalLinkage = "ward.D2",
                                  reps = 100,
                                  maxK = 10,
                                  seed = 123)

getCDFCurves(consClust)
getPAC(getCDFValues(consClust)$cdfDF)

#res <- M3C(mydata, des = desx, removeplots = TRUE, iters=25, objective='PAC', fsize=8, lthick=1, dotsize=1.25)
met_dat_4M3C <- l5_none_sepsis$sample_info
dat_4M3C <- t(l5_none_sepsis$stand)
#dimnames(dat_4M3C) <- lapply(dimnames(dat_4M3C), make.names)
colnames(met_dat_4M3C) <- gsub("sample", "ID", colnames(met_dat_4M3C))
met_dat_4M3C$ID <- make.names(met_dat_4M3C$ID)
colnames(dat_4M3C) <- gsub("_", ".", colnames(dat_4M3C), fixed = T)
colnames(dat_4M3C) <- paste0("samp.", colnames(dat_4M3C))
#dat_4M3C <- dat_4M3C[, 1:50]
#colnames(dat_4M3C) <- colnames(mydata_ori)

res <- M3C(dat_4M3C,
           removeplots = TRUE,
           iters=25,
           objective='PAC',
           fsize=8,
           lthick=1,
           dotsize=1.25)
cons_clust_ass <- data.frame(sample = gsub("samp.",
                                           "",
                                           colnames(dat_4M3C)),
                             group = factor(res$assignments))

write.csv(cons_clust_ass,
          file = sprintf("%scounts_l5_none_sepsis_unsup_subgroups.csv",
                         outDir))

comb_df$subgroup_cons_clust <- comb_df$group
comb_df$subgroup_cons_clust[comb_df$subgroup_cons_clust == "sepsis"] <- paste("sepsis",
                                                                              cons_clust_ass$group[match(comb_df$patient_ID[comb_df$subgroup_cons_clust == "sepsis"],
                                                                                                         gsub("\\_.*", "", cons_clust_ass$sample))],
                                                                              sep = "_")
table(comb_df$regrouping[comb_df$subgroup_cons_clust == "sepsis_1"])
table(comb_df$regrouping)
table(comb_df$subgroup)

cons_clust_ass$sample <- gsub(".", "_", cons_clust_ass$sample, fixed = T)

l5_none_sepsis_pca <- prcomp(l5_none_sepsis$stand, scale. = F, center = F)
l5_none_sepsis_OLs <- detectOLs_isfrst(l5_none_sepsis$stand,
                                       group_wise = T,
                                       samp_info = cons_clust_ass,
                                       groups_in = "group",
                                       quant_thrs = 0.99,
                                       pca_red_var = 0.95)

cons_clust_groups_pca <- plotUtils::plotPCA(l5_none_sepsis_pca,
                                            samp_info = cons_clust_ass,
                                            col = "group",
                                            fix_coord = F,
                                            labs = l5_none_sepsis_OLs$outliers)

ggsave(sprintf("%sl5_none_sepsis_cons_clust_groups_pca.pdf", outDir),
       plot = cons_clust_groups_pca,
       height = 4,
       width = 5)


comb_df$patient_ID
comb_df$sample <- comb_df$patient_ID
comb_df$sample <- paste(comb_df$sample, "none", sep = "_")

clinicl_subgroups_pca <- plotUtils::plotPCA(l5_none_sepsis_pca,
                                            samp_info = comb_df,
                                            col = "regrouping",
                                            fix_coord = F,
                                            labs = l5_none_sepsis_OLs$outliers)
ggsave(sprintf("%sl5_none_sepsis_clin_subgroups_pca.pdf", outDir),
       plot = clinicl_subgroups_pca,
       height = 4,
       width = 5)
