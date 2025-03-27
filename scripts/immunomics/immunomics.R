#remove.packages("plotUtils")
#rm(list = ls())
#detach("package:plotUtils")

if (!require("plotUtils", quietly = T)){
        devtools::install_github("guisantagui/plotUtils", upgrade = "never")
}
library(plotUtils)
library(ropls)

group_data_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/data/TIPS-Data_analysis/TIPS-Results_collection/TIPS_data_collection_2025-01-11/grouping_master_data_2025-01-10.RDS"
subgroup_data_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/results/immunomics/grouping_subgroups.csv"
#group_data_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/data/TIPS-Data_analysis/TIPS-Results_collection/TIPS_data_collection_2024-07-22/grouping_master_data.RDS"
immu_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/data/Immunomics_master_data_2025-01-10.RData"
outDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/results/immunomics/"
sepsis_subgroups_f <- sprintf("%s/counts_l5_none_sepsis_unsup_subgroups.csv",
                              outDir)
outDir <- add_slash_if_not(outDir)

create_dir_if_not(outDir)

load(immu_f)
group_data <- readRDS(group_data_f)
subgroup_data <- read.csv(subgroup_data_f, row.names = 1)
sepsis_subgroups <- read.csv(sepsis_subgroups_f, row.names = 1)
sepsis_subgroups$sample <- gsub("\\_.*", "", sepsis_subgroups$sample)
sepsis_subgroups$sample <- gsub("\\..*", "", sepsis_subgroups$sample)
group_data$subgroups <- group_data$group
group_data$subgroups[match(sepsis_subgroups$sample,
                           group_data$patient_ID)] <- paste(group_data$subgroups[match(sepsis_subgroups$sample,
                                                                                       group_data$patient_ID)],
                                                            sepsis_subgroups$group,
                                                            sep = "_")
group_data$clin_subgroup <- subgroup_data$regrouping[match(group_data$patient_ID,
                                                           subgroup_data$patient_ID)]
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
                immun_comm$subgroup <- patient_info$subgroup[match(immun_comm$patient_ID,
                                                                   patient_info$patient_ID)]
                immun_comm$clin_subgroup <- patient_info$clin_subgroup[match(immun_comm$patient_ID,
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

m <- matrix(1:20, nrow = 4)
dimnames(m) <- list(1:4, c("A", "B", "C", "D", "E"))
m[c(1, 3, 1), c(2, 4)] <- NA

na_mask <- is.na(m)

which(!na_mask, arr.ind = T)

m_stand_lst <- apply(m, 2, function(x) (x[!is.na(x)] - mean(x[!is.na(x)]))/sd(x[!is.na(x)]))

m_stand <- matrix(NA, nrow = nrow(m), ncol = ncol(m))
m_stand[!na_mask] <- unlist(m_stand_lst)

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
        slots <- slots[!grepl("OLs", slots)]
        for (s in slots){
                obj[[s]] <- obj[[s]][rownames(obj[[s]]) %in% obj$sample_info$sample, ]
        }
        return(obj)
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


# Plot distributions
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

lapply(imm_counts_list_parsed, function(x) plot(density(as.matrix(x$stand))))

plot(density(imm_counts_list_parsed_none$l1$stand$Bcell[imm_counts_list_parsed_none$l1$sample_info$group == "aid"]))
plot(density(imm_counts_list_parsed_none$l1$stand$Bcell[imm_counts_list_parsed_none$l1$sample_info$group == "sepsis"]))
plot(density(imm_counts_list_parsed_none$l1$stand$Bcell[imm_counts_list_parsed_none$l1$sample_info$group == "undefined"]))



lapply(imm_counts_list_parsed_none, function(x) plot(density(as.matrix(x$stand))))
lapply(imm_counts_list_parsed_LPS, function(x) plot(density(as.matrix(x$stand))))
lapply(imm_counts_list_parsed_R848, function(x) plot(density(as.matrix(x$stand))))

# PCAs
################################################################################
pca_counts_dir <- paste0(outDir, "pca/counts/")
create_dir_if_not(pca_counts_dir)

# Creates a list of multi-scoreplots from a list of PCAs
multPlot_pca_list <- function(pca_lst,
                              samp_info_list,
                              labs_list = rep(list(NULL), length(pca_lst)),
                              col = NULL,
                              nComps = 5,
                              biplot = T,
                              topNFeats = 5,
                              point_size = 3){
        plot_list <- mapply(function(pca, si, labs){
                plotUtils::doPCAMultiPlot(pca,
                                          samp_info = si,
                                          col = col,
                                          nComps = nComps,
                                          biplot = biplot,
                                          topNFeats = topNFeats,
                                          point_size = point_size,
                                          labs = labs)
                },
                pca_lst,
                samp_info_list,
                labs_list,
                SIMPLIFY = F)
        return(plot_list)
}

# Given a list of plots, saves them.
save_plotList <- function(plot_list, outDir, tag, height, width){
        mapply(function(p, n) save_plot(filename = sprintf("%s%s_%s.pdf",
                                                        outDir,
                                                        n,
                                                        tag),
                                     plot = p,
                                     height = height,
                                     width = width),
               plot_list,
               names(plot_list),
               SIMPLIFY = F)
}

# PCAs of all the treatments at once
imm_counts_pcas <- lapply(imm_counts_list_parsed,
                          function(x) prcomp(x$stand,
                                             scale. = F,
                                             center = F))

pca_counts_plts_group <- multPlot_pca_list(imm_counts_pcas,
                                    lapply(imm_counts_list_parsed,
                                           function(x) x$sample_info),
                                    col = "group",
                                    nComps = 5,
                                    point_size = .5)

pca_counts_plts_subgroup <- multPlot_pca_list(imm_counts_pcas,
                                              lapply(imm_counts_list_parsed,
                                                     function(x) x$sample_info),
                                              col = "subgroup",
                                              nComps = 5,
                                              point_size = .5)

pca_counts_plts_clinSubgroup <- multPlot_pca_list(imm_counts_pcas,
                                                  lapply(imm_counts_list_parsed,
                                                     function(x) x$sample_info),
                                                  col = "clin_subgroup",
                                                  nComps = 5,
                                                  point_size = .5)

save_plotList(pca_counts_plts_group,
              outDir = pca_counts_dir,
              tag = "counts_group_pca",
              height = 7,
              width = 7)

save_plotList(pca_counts_plts_subgroup,
              outDir = pca_counts_dir,
              tag = "counts_subgroup_pca",
              height = 7,
              width = 7)

save_plotList(pca_counts_plts_clinSubgroup,
              outDir = pca_counts_dir,
              tag = "counts_clin_subgroup_pca",
              height = 7,
              width = 7)

pca_counts_plts_treat <- multPlot_pca_list(imm_counts_pcas,
                                           lapply(imm_counts_list_parsed,
                                                  function(x) x$sample_info),
                                           col = "treatment",
                                           nComps = 5,
                                           point_size = .5)

save_plotList(pca_counts_plts_treat,
              outDir = pca_counts_dir,
              tag = "counts_treat_pca",
              height = 7,
              width = 7)

# PCAs of each treatment separately.

# No treatment
pca_counts_none_dir <- paste0(pca_counts_dir, "none/")
create_dir_if_not(pca_counts_none_dir)

imm_counts_pcas_none <- lapply(imm_counts_list_parsed_none,
                        function(x) prcomp(x$stand,
                                           scale. = F,
                                           center = F))

pca_counts_plts_none <- multPlot_pca_list(imm_counts_pcas_none,
                                          lapply(imm_counts_list_parsed_none,
                                                 function(x) x$sample_info),
                                   col = "group",
                                   nComps = 5,
                                   point_size = .5)
save_plotList(pca_counts_plts_none,
              outDir = pca_counts_none_dir,
              tag = "counts_none_mult_pca",
              height = 7,
              width = 7)

l5_counts_none_pca <- plotPCA(imm_counts_pcas_none$l5,
                              samp_info = imm_counts_list_parsed_none$l5$sample_info,
                              col = "group",
                              fix_coord = F,
                              biplot = T,
                              topNFeats = 5)

save_plot(filename = sprintf("%sl5_counts_none_pca.pdf",
                             pca_counts_none_dir),
          plot = l5_counts_none_pca,
          height = 4,
          width = 5)

l5_counts_none_pca_subg <- plotPCA(imm_counts_pcas_none$l5,
                                   samp_info = imm_counts_list_parsed_none$l5$sample_info,
                                   col = "subgroup",
                                   fix_coord = F,
                                   biplot = T,
                                   topNFeats = 5)

save_plot(filename = sprintf("%sl5_counts_none_pca_subgroup.pdf",
                             pca_counts_none_dir),
          plot = l5_counts_none_pca_subg,
          height = 4,
          width = 5)

l5_counts_none_pca_clinSubg <- plotPCA(imm_counts_pcas_none$l5,
                                       samp_info = imm_counts_list_parsed_none$l5$sample_info,
                                       col = "clin_subgroup",
                                       fix_coord = F,
                                       biplot = T,
                                       topNFeats = 5)

save_plot(filename = sprintf("%sl5_counts_none_pca_clin_subgroup.pdf",
                             pca_counts_none_dir),
          plot = l5_counts_none_pca_clinSubg,
          height = 4,
          width = 5)

# LPS
pca_counts_LPS_dir <- paste0(pca_counts_dir, "LPS/")
create_dir_if_not(pca_counts_LPS_dir)

imm_counts_pcas_LPS <- lapply(imm_counts_list_parsed_LPS,
                       function(x) prcomp(x$stand,
                                          scale. = F,
                                          center = F))

pca_counts_plts_LPS <- multPlot_pca_list(imm_counts_pcas_LPS,
                                         lapply(imm_counts_list_parsed_LPS,
                                                function(x) x$sample_info),
                                  col = "group",
                                  nComps = 5,
                                  point_size = .5)
save_plotList(pca_counts_plts_LPS,
              outDir = pca_counts_LPS_dir,
              tag = "LPS_mult_pca",
              height = 7,
              width = 7)

l5_counts_LPS_pca <- plotPCA(imm_counts_pcas_LPS$l5,
                             samp_info = imm_counts_list_parsed_LPS$l5$sample_info,
                             col = "group",
                             fix_coord = F,
                             biplot = T,
                             topNFeats = 5)

save_plot(filename = sprintf("%sl5_counts_LPS_pca.pdf",
                             pca_counts_LPS_dir),
          plot = l5_counts_LPS_pca,
          height = 4,
          width = 5)

l5_counts_LPS_pca_subg <- plotPCA(imm_counts_pcas_LPS$l5,
                                  samp_info = imm_counts_list_parsed_LPS$l5$sample_info,
                                  col = "subgroup",
                                  fix_coord = F,
                                  biplot = T,
                                  topNFeats = 5)

save_plot(filename = sprintf("%sl5_LPS_pca_subgroup.pdf",
                             pca_counts_LPS_dir),
          plot = l5_counts_LPS_pca_subg,
          height = 4,
          width = 5)

l5_counts_LPS_pca_clinSubg <- plotPCA(imm_counts_pcas_LPS$l5,
                                      samp_info = imm_counts_list_parsed_LPS$l5$sample_info,
                                      col = "clin_subgroup",
                                      fix_coord = F,
                                      biplot = T,
                                      topNFeats = 5)

save_plot(filename = sprintf("%sl5_counts_LPS_pca_clin_subgroup.pdf",
                             pca_counts_LPS_dir),
          plot = l5_counts_LPS_pca_clinSubg,
          height = 4,
          width = 5)

# R848
pca_counts_R848_dir <- paste0(pca_counts_dir, "R848/")
create_dir_if_not(pca_counts_R848_dir)

imm_counts_pcas_R848 <- lapply(imm_counts_list_parsed_R848,
                        function(x) prcomp(x$stand,
                                           scale. = F,
                                           center = F))

pca_counts_plts_R848 <- multPlot_pca_list(imm_counts_pcas_R848,
                                          lapply(imm_counts_list_parsed_R848,
                                         function(x) x$sample_info),
                                  col = "group",
                                  nComps = 5,
                                  point_size = .5)
save_plotList(pca_counts_plts_R848,
              outDir = pca_counts_R848_dir,
              tag = "R848_mult_pca",
              height = 7,
              width = 7)

l5_counts_R848_pca <- plotPCA(imm_counts_pcas_R848$l5,
                             samp_info = imm_counts_list_parsed_R848$l5$sample_info,
                             col = "group",
                             fix_coord = F,
                             biplot = T,
                             topNFeats = 5)

save_plot(filename = sprintf("%sl5_counts_R848_pca.pdf",
                             pca_counts_R848_dir),
          plot = l5_counts_R848_pca,
          height = 4,
          width = 5)

l5_counts_R848_pca_subg <- plotPCA(imm_counts_pcas_R848$l5,
                                   samp_info = imm_counts_list_parsed_R848$l5$sample_info,
                                   col = "subgroup",
                                   fix_coord = F,
                                   biplot = T,
                                   topNFeats = 5)

save_plot(filename = sprintf("%sl5_R848_pca_subgroup.pdf",
                             pca_counts_R848_dir),
          plot = l5_counts_R848_pca_subg,
          height = 4,
          width = 5)

l5_counts_R848_pca_clinSubg <- plotPCA(imm_counts_pcas_R848$l5,
                                       samp_info = imm_counts_list_parsed_R848$l5$sample_info,
                                       col = "clin_subgroup",
                                       fix_coord = F,
                                       biplot = T,
                                       topNFeats = 5)

save_plot(filename = sprintf("%sl5_counts_R848_pca_clin_subgroup.pdf",
                             pca_counts_R848_dir),
          plot = l5_counts_R848_pca_clinSubg,
          height = 4,
          width = 5)

# OPLS-DAs
################################################################################

do_opls_groups_vs_rest <- function(obj,
                                   slot = "stand",
                                   orthoI = 2,
                                   permI = 100,
                                   plot_dir = NULL,
                                   tag = NULL){
        #obj <- imm_states_list_parsed_LPS$l1
        #slot <- "stand"
        #orthoI = 3
        #permI = 100
        #plot_dir = paste0(outDir, "opls/states/LPS/plots/")
        #tag = "LPS"
        group_vec <- obj$sample_info$group[match(rownames(obj[[slot]]),
                                                 obj$sample_info$sample)]
        groups <- unique(group_vec)
        opls_list <- list()
        if (slot == "stand"){
                scaling <- "none"
        }else{
                scaling <- "standard"
        }
        opls_fun <- function(dat,
                             y,
                             orthoI,
                             permI,
                             scaling,
                             plot_dir,
                             list_name,
                             tag){
                if (!is.null(plot_dir)){
                        if (!dir.exists(plot_dir)){
                                warning(sprintf("%s doesn't exist. Creating it...",
                                                plot_dir),
                                        call. = F)
                                dir.create(plot_dir, recursive = T)
                        }
                        if (!is.null(tag)){
                                out_file <- sprintf("%sopls_%s_%s.pdf",
                                                    plot_dir,
                                                    tag,
                                                    list_name)
                        }else{
                                out_file <- sprintf("%sopls_%s.pdf",
                                                    plot_dir,
                                                    list_name)
                        }
                        pdf(file = out_file, height = 5, width = 5)
                        opls_da <- opls(dat,
                                        y = y,
                                        orthoI = orthoI,
                                        predI = 1,
                                        permI = permI,
                                        scaleC = scaling)
                        dev.off()
                }else{
                        opls_da <- opls(dat,
                                        y = y,
                                        orthoI = orthoI,
                                        predI = 1,
                                        permI = permI,
                                        scaleC = scaling)
                }
                return(opls_da)
        }
        if (length(groups) == 2){
                y <- factor(group_vec)
                list_name <- sprintf("%s_vs_%s",
                                     levels(y)[1],
                                     levels(y)[2])
                opls_da <- opls_fun(obj[[slot]],
                                    y = y,
                                    orthoI = orthoI,
                                    permI = permI,
                                    scaling = scaling,
                                    plot_dir = plot_dir,
                                    list_name = list_name,
                                    tag = tag)
                opls_list[[list_name]] <- opls_da
        }else{
                for (g in groups){
                        y <- factor(c("rest", g)[as.numeric(group_vec == g) + 1],
                                    levels = c("rest", g))
                        list_name <- sprintf("%s_vs_rest", g)
                        opls_da <- opls_fun(obj[[slot]],
                                            y = y,
                                            orthoI = orthoI,
                                            permI = permI,
                                            scaling = scaling,
                                            plot_dir = plot_dir,
                                            list_name,
                                            tag = tag)
                        opls_list[[list_name]] <- opls_da
                }
        }
        return (opls_list)
}







do_opls_list <- function(obj_list, slot, orthoI, permI, plot_dir = NULL, tag){
        opls_res <- mapply(function(x, n){
                do_opls_groups_vs_rest(x,
                                       slot = slot,
                                       orthoI = orthoI,
                                       permI = permI,
                                       plot_dir = plot_dir,
                                       tag = sprintf("%s_%s", tag, n))
        },
        obj_list,
        names(obj_list),
        SIMPLIFY = F)
        return(opls_res)
}

get_opls_sign <- function(obj, vip_thrshld = 0, vip_sort = T, save_dir = NULL){
        list_of_lists <- list()
        for (i in seq_along(obj)){
                n <- names(obj)[i]
                l <- obj[[i]]
                comps <- names(l)
                df_list <- list()
                for (comp in comps){
                        out_df <- data.frame(variable = names(l[[comp]]@vipVn),
                                             VIP = l[[comp]]@vipVn,
                                             loading = l[[comp]]@loadingMN[, 1])
                        out_df <- out_df[out_df$VIP >= vip_thrshld, ]
                        if (vip_sort){
                                out_df <- out_df[order(out_df$VIP, decreasing = T), ]
                        }
                        df_list[[comp]] <- out_df
                        if (!is.null(save_dir)){
                                if (!dir.exists(save_dir)){
                                        warning(sprintf("%s doesn't exist. Creating it...",
                                                        save_dir),
                                                call. = F)
                                        dir.create(save_dir, recursive = T)
                                }
                                write.csv(out_df,
                                          file = sprintf("%s%s_%s.csv", save_dir, n, comp))
                        }
                }
                list_of_lists[[n]] <- df_list
        }
        return(list_of_lists)
}

get_summary_df <- function(opls_list){
        #opls_list <- opls_states_LPS_aid_vs_sepsis
        df_out <- data.frame(matrix(nrow = 0, ncol = 10,
                                    dimnames = list(NULL,
                                                    c("level",
                                                      "comparison",
                                                      "R2X(cum)",
                                                      "R2Y(cum)",
                                                      "Q2(cum)",
                                                      "RMSEE",
                                                      "pre",
                                                      "ort",
                                                      "pR2Y",
                                                      "pQ2"))))
        for (i in seq_along(opls_list)){
                l <- names(opls_list)[i]
                level <- opls_list[[l]]
                for (comp in names(level)){
                        comp <- names(level)[1]
                        df <- level[[comp]]@summaryDF
                        df$comparison <- comp
                        df$level <- l
                        df <- df[, c(10, 9, 1:8)]
                        df_out <- rbind.data.frame(df_out, df)
                }
        }
        rownames(df_out) <- 1:nrow(df_out)
        return(df_out)
}

opls_counts_LPS <- do_opls_list(imm_counts_list_parsed_LPS,
                                slot = "stand",
                                orthoI = 3,
                                permI = 100,
                                plot_dir = paste0(outDir, "opls/counts/LPS/plots/"),
                                tag = "LPS")

opls_counts_R848 <- do_opls_list(imm_counts_list_parsed_R848,
                                 slot = "stand",
                                 orthoI = 3,
                                 permI = 100,
                                 plot_dir = paste0(outDir, "opls/counts/R848/plots/"),
                                 tag = "R848")

opls_counts_none <- do_opls_list(imm_counts_list_parsed_none,
                                 slot = "stand",
                                 orthoI = 3,
                                 permI = 100,
                                 plot_dir = paste0(outDir, "opls/counts/none/plots/"),
                                 tag = "none")



opls_counts_LPS_sign <- get_opls_sign(opls_counts_LPS, vip_thrshld = 1,
                                      save_dir = paste0(outDir,
                                                        "opls/counts/LPS/csv/"))
opls_counts_R848_sign <- get_opls_sign(opls_counts_R848, vip_thrshld = 1,
                                       save_dir = paste0(outDir,
                                                         "opls/counts/R848/csv/"))
opls_counts_none_sign <- get_opls_sign(opls_counts_none, vip_thrshld = 1,
                                       save_dir = paste0(outDir,
                                                         "opls/counts/none/csv/"))
# Do OPLS of aid vs sepsis
imm_counts_list_parsed_LPS_no_undef <- lapply(imm_counts_list_parsed_LPS,
                                              filt_samps,
                                              keep_var = "group",
                                              keep_val = c("aid", "sepsis"))
imm_counts_list_parsed_R848_no_undef <- lapply(imm_counts_list_parsed_R848,
                                               filt_samps,
                                               keep_var = "group",
                                               keep_val = c("aid", "sepsis"))
imm_counts_list_parsed_none_no_undef <- lapply(imm_counts_list_parsed_none,
                                               filt_samps,
                                               keep_var = "group",
                                               keep_val = c("aid", "sepsis"))


opls_counts_LPS_aid_vs_sepsis <- do_opls_list(imm_counts_list_parsed_LPS_no_undef,
                                              slot = "stand",
                                              orthoI = 3,
                                              permI = 100,
                                              plot_dir = paste0(outDir,
                                                                "opls/counts/LPS/plots/"),
                                              tag = "LPS")

opls_counts_R848_aid_vs_sepsis <- do_opls_list(imm_counts_list_parsed_R848_no_undef,
                                               slot = "stand",
                                               orthoI = 3,
                                               permI = 100,
                                               plot_dir = paste0(outDir,
                                                                 "opls/counts/R848/plots/"),
                                               tag = "R848")

opls_counts_none_aid_vs_sepsis <- do_opls_list(imm_counts_list_parsed_none_no_undef,
                                               slot = "stand",
                                               orthoI = 3,
                                               permI = 100,
                                               plot_dir = paste0(outDir, "opls/counts/none/plots/"),
                                               tag = "none")

opls_summ_cnts_LPS_aid_vs_sep <- get_summary_df(opls_counts_LPS_aid_vs_sepsis)
opls_summ_cnts_R848_aid_vs_sep <- get_summary_df(opls_counts_R848_aid_vs_sepsis)
opls_summ_cnts_none_aid_vs_sep <- get_summary_df(opls_counts_none_aid_vs_sepsis)

opls_summ_cnts_LPS_aid_vs_sep$treatment <- "LPS"
opls_summ_cnts_R848_aid_vs_sep$treatment <- "R848"
opls_summ_cnts_none_aid_vs_sep$treatment <- "none"

opls_summ_cnts_aid_vs_sep <- rbind.data.frame(opls_summ_cnts_none_aid_vs_sep,
                                              opls_summ_cnts_LPS_aid_vs_sep,
                                              opls_summ_cnts_R848_aid_vs_sep)
opls_summ_cnts_aid_vs_sep <- opls_summ_cnts_aid_vs_sep[, c(1:6, 9, 10, 11)]
opls_summ_cnts_aid_vs_sep[, c(ncol(opls_summ_cnts_aid_vs_sep),
                              1:(ncol(opls_summ_cnts_aid_vs_sep) - 1))]

################################################################################
################################################################################
# States analysis                                                              #
################################################################################
################################################################################

imm_states_list <- list(l1 = immunomics_states_level1,
                        l2 = immunomics_states_level2,
                        l3 = immunomics_states_level3,
                        l4 = immunomics_states_level4,
                        l5 = immunomics_states_level5)

# Parse the files
imm_states_list_parsed <- lapply(imm_states_list,
                                 parse_immun,
                                 patient_info = group_data,
                                 type = "states")

# Filter to keep samples we have
imm_states_list_parsed <- lapply(imm_states_list_parsed,
                                 function(x) filt_samps(x,
                                                        keep_samps = x$sample_info$sample[!is.na(x$sample_info$group)]))

# Functions for outlier detection

# Standard deviation method
detectOLs_SD <- function(citoNum, sdThrshld = 2, groupWise = F){
        if(!require(stringr)) install.packages("stringr")
        library(stringr)
        if(groupWise){
                unCases <- unique(str_extract(rownames(citoNum), "[^_]*_[^_]*"))
                logicDF <- data.frame(matrix(nrow = 0,
                                             ncol = ncol(citoNum),
                                             dimnames = list(NULL,
                                                             colnames(citoNum))))
                for(c in unCases){
                        citCase <- citoNum[grep(c, rownames(citoNum)), ]
                        citCase <- apply(citCase, 2, function(x) (x - mean(x))/sd(x))
                        logCase <- citCase > sdThrshld | citCase < -sdThrshld
                        logicDF <- rbind.data.frame(logicDF, logCase)
                }
        }else{
                scaledCentered <- apply(citoNum, 2, function(x) (x - mean(x))/sd(x))
                logicDF <- scaledCentered > sdThrshld | scaledCentered < -sdThrshld
        }
        OLs <- apply(logicDF, 1, any)
        OLs <- names(OLs)[OLs]
        return(OLs)
}

if (!require("isotree", quietly = T)){
        install.packages("isotree")
}
library(isotree)

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

# IQR method
detectOLs_IQR <- function(m, k = 1.5, groupWise = F, logDF = F){
        m <- imm_states_list_parsed$l1$stand
        if(!require(stringr)) install.packages("stringr")
        library(stringr)
        getLogic <- function(m, k){
                logicFun <- data.frame(matrix(nrow = nrow(m),
                                              ncol = 0,
                                              dimnames = list(rownames(m),
                                                              NULL)))
                for(i in 1:ncol(m)){
                        name <- colnames(m)[i]
                        quants <- quantile(m[, i], probs = c(.25, .75))
                        iqr <- quants[1] - quants[2]
                        cutoff <- iqr * k
                        interval <- c(quants[2] - cutoff,
                                      quants[1] + cutoff)
                        log2Bind <- data.frame(matrix(m[, i] < interval[1] | m[, i] > interval[2],
                                                      nrow = nrow(m),
                                                      ncol = 1,
                                                      dimnames = list(rownames(m),
                                                                      name)))
                        logicFun <- cbind.data.frame(logicFun, log2Bind)
                }
                return(logicFun)
        }
        
        if(groupWise){
                logicDF <- data.frame(matrix(nrow = 0,
                                             ncol = ncol(m),
                                             dimnames = list(NULL,
                                                             colnames(m))))
                unCases <- unique(str_extract(rownames(m), "[^_]*_[^_]*"))
                for(c in unCases){
                        citCase <- m[grep(c, rownames(m)), ]
                        subLog <- getLogic(citCase, k)
                        logicDF <- rbind.data.frame(logicDF,
                                                    subLog)
                }
        }else{
                logicDF <- getLogic(m, k)
        }
        if(logDF){
                print(logicDF)
        }
        OLs <- apply(logicDF, 1, any)
        OLs <- names(OLs)[OLs]
        return(OLs)
}



library(impute)

impute_list <- function(obj, slot,
                        max_in_row,
                        max_in_col,
                        k = 10,
                        filt_from_samp_info = T){
        for (i in seq_along(obj)){
                l <- obj[[i]]
                keep_col <- apply(l[[slot]],
                                  2,
                                  function(x) sum(is.na(x))/length(x) <= max_in_col)
                keep_row <- apply(l[[slot]],
                                  1,
                                  function(x) sum(is.na(x))/length(x) <= max_in_row)
                l[[slot]] <- l[[slot]][keep_row, keep_col]
                imp_flags <- is.na(as.matrix(l[[slot]]))
                l[[slot]] <- as.data.frame(impute.knn(as.matrix(l[[slot]]),
                                                      rowmax = 1,
                                                      colmax = 1,
                                                      k = k)$data)
                obj[[i]][[sprintf("%s_imp", slot)]] <- l[[slot]]
                obj[[i]][[sprintf("%s_imp_flags", slot)]] <- imp_flags
                if (filt_from_samp_info){
                        obj[[i]]$sample_info <- obj[[i]]$sample_info[keep_row, ]
                }
        }
        return(obj)
}

# Applies log2 over a list
do_log2_list <- function(obj, slot, add_base = 0){
        if (!is.null(add_base)){
                add_base <- rep(add_base, length(obj))
        }else{
                add_base <- sapply(obj,
                                   function(x) min(x[[slot]][x[[slot]] != 0 & !is.na(x[[slot]])]))
        }
        obj <- mapply(function(x, m){
                do_log2(x,
                        slot = slot,
                        add_base = m)
        },
        obj,
        add_base,
        SIMPLIFY = F)
        return(obj)
}

imm_states_list_parsed <- impute_list(imm_states_list_parsed,
                                      slot = "data",
                                      max_in_row = 0.2,
                                      max_in_col = 0.8)

imm_states_list_parsed <- do_log2_list(imm_states_list_parsed,
                                       slot = "data_imp",
                                       add_base = NULL)

imm_states_list_parsed <- lapply(imm_states_list_parsed,
                                 do_stand,
                                 slot = "log2")

plot(density(as.matrix(imm_states_list_parsed$l1$stand)))

plot(density(as.matrix(log2(imm_states_list_parsed$l1$data_imp + min(imm_states_list_parsed$l1$data_imp[imm_states_list_parsed$l1$data_imp != 0])))))

detectOLs_isfrst(imm_states_list_parsed$l1$stand, group_wise = T, samp_info = imm_states_list_parsed$l1$sample_info,
                 groups_in = "group")


# PCAs
################################################################################
pca_states_dir <- paste0(outDir, "pca/states/")
create_dir_if_not(pca_states_dir)

# Apply isolation forest outlier detection over list object
detectOLs_isfrst_list <- function(obj,
                                  slot,
                                  quant_thrs = 0.95,
                                  group_wise = F,
                                  groups_in,
                                  pca_red_var){
        for (i in seq_along(obj)){
                ols <- detectOLs_isfrst(obj[[i]][[slot]],
                                        quant_thrs = quant_thrs,
                                        group_wise = group_wise,
                                        samp_info = obj[[i]]$sample_info,
                                        groups_in = groups_in,
                                        pca_red_var = pca_red_var)
                obj[[i]][[sprintf("%s_OLs", slot)]] <- ols
        }
        return(obj)
}

# Detect outliers
imm_states_list_parsed <- detectOLs_isfrst_list(imm_states_list_parsed,
                                                slot = "stand",
                                                quant_thrs = 0.95,
                                                group_wise = T,
                                                groups_in = "group",
                                                pca_red_var = 0.95)
imm_states_list_parsed$l5$stand_OLs$anomaly_scores
imm_states_list_parsed$l5$stand_OLs$outliers
# PCAs of all the treatments at once
imm_states_pcas <- lapply(imm_states_list_parsed,
                          function(x) prcomp(x$stand,
                                             scale. = F,
                                             center = F))

pca_imp <- summary(imm_states_pcas$l5)$importance


pca_states_plts_group <- multPlot_pca_list(imm_states_pcas,
                                           lapply(imm_states_list_parsed,
                                                  function(x) x$sample_info),
                                           col = "group",
                                           nComps = 5,
                                           point_size = 0.5,
                                           labs_list = lapply(imm_states_list_parsed,
                                                              function(x) x$stand_OLs$outliers))

pca_states_plts_clinSubgroup <- multPlot_pca_list(imm_states_pcas,
                                                  lapply(imm_states_list_parsed,
                                                         function(x) x$sample_info),
                                                  col = "clin_subgroup",
                                                  nComps = 5,
                                                  point_size = .5,
                                                  labs_list = lapply(imm_states_list_parsed,
                                                                     function(x) x$stand_OLs$outliers))

save_plotList(pca_states_plts_group,
              outDir = pca_states_dir,
              tag = "states_group_pca",
              height = 7,
              width = 7)

save_plotList(pca_states_plts_clinSubgroup,
              outDir = pca_states_dir,
              tag = "states_clin_subgroup_pca",
              height = 7,
              width = 7)


pca_states_plts_treat <- multPlot_pca_list(imm_states_pcas,
                                           lapply(imm_states_list_parsed,
                                                  function(x) x$sample_info),
                                           col = "treatment",
                                           nComps = 5,
                                           point_size = 0.5,
                                           labs_list = lapply(imm_states_list_parsed,
                                                              function(x) x$stand_OLs$outliers))

save_plotList(pca_states_plts_treat,
              outDir = pca_states_dir,
              tag = "states_treat_pca",
              height = 7,
              width = 7)

# PCAs of each treatment separately.

# Plot distributions
imm_states_list_parsed_none <- lapply(imm_states_list_parsed,
                                      function(x) filt_samps(x,
                                                             keep_var = "treatment",
                                                             keep_val = "none"))
imm_states_list_parsed_LPS <- lapply(imm_states_list_parsed,
                                     function(x) filt_samps(x,
                                                            keep_var = "treatment",
                                                            keep_val = "LPS"))
imm_states_list_parsed_R848 <- lapply(imm_states_list_parsed,
                                      function(x) filt_samps(x,
                                                             keep_var = "treatment",
                                                             keep_val = "R848"))

# Stand again
imm_states_list_parsed_none <- lapply(imm_states_list_parsed_none,
                                      do_stand,
                                      slot = "log2")
imm_states_list_parsed_LPS <- lapply(imm_states_list_parsed_LPS,
                                     do_stand,
                                     slot = "log2")
imm_states_list_parsed_R848 <- lapply(imm_states_list_parsed_R848,
                                      do_stand,
                                      slot = "log2")

# Find outliers for each dataset, separated
imm_states_list_parsed_none <- detectOLs_isfrst_list(imm_states_list_parsed_none,
                                                     slot = "stand",
                                                     quant_thrs = 0.95,
                                                     group_wise = T,
                                                     groups_in = "group",
                                                     pca_red_var = 0.95)

imm_states_list_parsed_LPS <- detectOLs_isfrst_list(imm_states_list_parsed_LPS,
                                                     slot = "stand",
                                                     quant_thrs = 0.95,
                                                     group_wise = T,
                                                     groups_in = "group",
                                                     pca_red_var = 0.95)

imm_states_list_parsed_R848 <- detectOLs_isfrst_list(imm_states_list_parsed_R848,
                                                     slot = "stand",
                                                     quant_thrs = 0.95,
                                                     group_wise = T,
                                                     groups_in = "group",
                                                     pca_red_var = 0.95)

# No treatment
pca_states_none_dir <- paste0(pca_states_dir, "none/")
create_dir_if_not(pca_states_none_dir)

imm_states_pcas_none <- lapply(imm_states_list_parsed_none,
                               function(x) prcomp(x$stand,
                                                  scale. = F,
                                                  center = F))

pca_states_plts_none <- multPlot_pca_list(imm_states_pcas_none,
                                          lapply(imm_states_list_parsed_none,
                                                 function(x) x$sample_info),
                                          col = "group",
                                          nComps = 5,
                                          point_size = 0.5,
                                          labs_list = lapply(imm_states_list_parsed_none,
                                                             function(x) x$stand_OLs$outliers))
save_plotList(pca_states_plts_none,
              outDir = pca_states_none_dir,
              tag = "states_none_mult_pca",
              height = 7,
              width = 7)


pca_states_plts_none_clinSubgrp <- multPlot_pca_list(imm_states_pcas_none,
                                                     lapply(imm_states_list_parsed_none,
                                                            function(x) x$sample_info),
                                                     col = "clin_subgroup",
                                                     nComps = 5,
                                                     point_size = 0.5,
                                                     labs_list = lapply(imm_states_list_parsed_none,
                                                                        function(x) x$stand_OLs$outliers))
save_plotList(pca_states_plts_none_clinSubgrp,
              outDir = pca_states_none_dir,
              tag = "states_none_mult_clinSubgrp_pca",
              height = 7,
              width = 7)

l5_states_none_pca <- plotPCA(imm_states_pcas_none$l5,
                              samp_info = imm_states_list_parsed_none$l5$sample_info,
                              col = "group",
                              fix_coord = F,
                              biplot = T,
                              topNFeats = 5,
                              labs = imm_states_list_parsed_none$l5$stand_OLs$outliers)

save_plot(filename = sprintf("%sl5_states_none_pca.pdf",
                             pca_states_none_dir),
          plot = l5_states_none_pca,
          height = 4,
          width = 5)

# LPS
pca_states_LPS_dir <- paste0(pca_states_dir, "LPS/")
create_dir_if_not(pca_states_LPS_dir)

imm_states_pcas_LPS <- lapply(imm_states_list_parsed_LPS,
                              function(x) prcomp(x$stand,
                                                 scale. = F,
                                                 center = F))

pca_states_plts_LPS <- multPlot_pca_list(imm_states_pcas_LPS,
                                         lapply(imm_states_list_parsed_LPS,
                                                function(x) x$sample_info),
                                         col = "group",
                                         nComps = 5,
                                         point_size = .5,
                                         labs_list = lapply(imm_states_list_parsed_LPS,
                                                            function(x) x$stand_OLs$outliers))
save_plotList(pca_states_plts_LPS,
              outDir = pca_states_LPS_dir,
              tag = "LPS_mult_pca",
              height = 7,
              width = 7)

pca_states_plts_LPS_clinSubgrp <- multPlot_pca_list(imm_states_pcas_LPS,
                                                    lapply(imm_states_list_parsed_LPS,
                                                           function(x) x$sample_info),
                                                    col = "clin_subgroup",
                                                    nComps = 5,
                                                    point_size = 0.5,
                                                    labs_list = lapply(imm_states_list_parsed_LPS,
                                                                       function(x) x$stand_OLs$outliers))
save_plotList(pca_states_plts_LPS_clinSubgrp,
              outDir = pca_states_LPS_dir,
              tag = "states_LPS_mult_clinSubgrp_pca",
              height = 7,
              width = 7)

l5_states_LPS_pca <- plotPCA(imm_states_pcas_LPS$l5,
                              samp_info = imm_states_list_parsed_LPS$l5$sample_info,
                              col = "group",
                              fix_coord = F,
                              biplot = T,
                              topNFeats = 5,
                              labs = imm_states_list_parsed_LPS$l5$stand_OLs$outliers)

save_plot(filename = sprintf("%sl5_states_LPS_pca.pdf",
                             pca_states_LPS_dir),
          plot = l5_states_LPS_pca,
          height = 4,
          width = 5)

# R848
pca_states_R848_dir <- paste0(pca_states_dir, "R848/")
create_dir_if_not(pca_states_R848_dir)

imm_states_pcas_R848 <- lapply(imm_states_list_parsed_R848,
                               function(x) prcomp(x$stand,
                                                  scale. = F,
                                                  center = F))

pca_states_plts_R848 <- multPlot_pca_list(imm_states_pcas_R848,
                                          lapply(imm_states_list_parsed_R848,
                                                 function(x) x$sample_info),
                                          col = "group",
                                          nComps = 5,
                                          point_size = .5,
                                          labs_list = lapply(imm_states_list_parsed_R848,
                                                             function(x) x$stand_OLs$outliers))
save_plotList(pca_states_plts_R848,
              outDir = pca_states_R848_dir,
              tag = "R848_mult_pca",
              height = 7,
              width = 7)

pca_states_plts_R848_clinSubgrp <- multPlot_pca_list(imm_states_pcas_R848,
                                                     lapply(imm_states_list_parsed_R848,
                                                            function(x) x$sample_info),
                                                     col = "clin_subgroup",
                                                     nComps = 5,
                                                     point_size = 0.5,
                                                     labs_list = lapply(imm_states_list_parsed_R848,
                                                                        function(x) x$stand_OLs$outliers))
save_plotList(pca_states_plts_R848_clinSubgrp,
              outDir = pca_states_R848_dir,
              tag = "states_R848_mult_clinSubgrp_pca",
              height = 7,
              width = 7)


l5_states_R848_pca <- plotPCA(imm_states_pcas_R848$l5,
                              samp_info = imm_states_list_parsed_R848$l5$sample_info,
                              col = "group",
                              fix_coord = F,
                              biplot = T,
                              topNFeats = 5,
                              labs = imm_states_list_parsed_R848$l5$stand_OLs$outliers)

save_plot(filename = sprintf("%sl5_states_R848_pca.pdf",
                             pca_states_R848_dir),
          plot = l5_states_R848_pca,
          height = 4,
          width = 5)

# Filter to keep only sepsis and see if we can see clustering according to 
# clinical subgrouping
imm_states_list_parsed_none_sepsis <- lapply(imm_states_list_parsed_none,
                                             function(x) filt_samps(x,
                                                                    keep_var = "group",
                                                                    keep_val = "sepsis"))
imm_states_list_parsed_LPS_sepsis <- lapply(imm_states_list_parsed_LPS,
                                            function(x) filt_samps(x,
                                                                   keep_var = "group",
                                                                   keep_val = "sepsis"))
imm_states_list_parsed_R848_sepsis <- lapply(imm_states_list_parsed_R848,
                                             function(x) filt_samps(x,
                                                                    keep_var = "group",
                                                                    keep_val = "sepsis"))

# Stand again
imm_states_list_parsed_none_sepsis <- lapply(imm_states_list_parsed_none_sepsis,
                                             do_stand,
                                             slot = "log2")
imm_states_list_parsed_LPS_sepsis <- lapply(imm_states_list_parsed_LPS_sepsis,
                                            do_stand,
                                            slot = "log2")
imm_states_list_parsed_R848_sepsis <- lapply(imm_states_list_parsed_R848_sepsis,
                                             do_stand,
                                             slot = "log2")

imm_states_pcas_none_sepsis <- lapply(imm_states_list_parsed_none_sepsis,
                                      function(x) prcomp(x$stand,
                                                         scale. = F,
                                                         center = F))

imm_states_pcas_LPS_sepsis <- lapply(imm_states_list_parsed_LPS_sepsis,
                                     function(x) prcomp(x$stand,
                                                        scale. = F,
                                                        center = F))

imm_states_pcas_R848_sepsis <- lapply(imm_states_list_parsed_R848_sepsis,
                                      function(x) prcomp(x$stand,
                                                         scale. = F,
                                                         center = F))

pca_states_plts_none_sepsis <- multPlot_pca_list(imm_states_pcas_none_sepsis,
                                                 lapply(imm_states_list_parsed_none_sepsis,
                                                        function(x) x$sample_info),
                                                 col = "clin_subgroup",
                                                 nComps = 5,
                                                 point_size = .5,
                                                 biplot = F)

pca_states_plts_LPS_sepsis <- multPlot_pca_list(imm_states_pcas_LPS_sepsis,
                                                lapply(imm_states_list_parsed_LPS_sepsis,
                                                       function(x) x$sample_info),
                                                col = "clin_subgroup",
                                                nComps = 5,
                                                point_size = .5,
                                                biplot = F)

pca_states_plts_R848_sepsis <- multPlot_pca_list(imm_states_pcas_R848_sepsis,
                                                 lapply(imm_states_list_parsed_R848_sepsis,
                                                        function(x) x$sample_info),
                                                 col = "clin_subgroup",
                                                 nComps = 5,
                                                 point_size = .5,
                                                 biplot = F)

save_plotList(pca_states_plts_none_sepsis,
              outDir = pca_states_none_dir,
              tag = "states_none_mult_sepsis_clinSubgrp_pca",
              height = 7,
              width = 7)

save_plotList(pca_states_plts_LPS_sepsis,
              outDir = pca_states_LPS_dir,
              tag = "states_LPS_mult_sepsis_clinSubgrp_pca",
              height = 7,
              width = 7)

save_plotList(pca_states_plts_R848_sepsis,
              outDir = pca_states_R848_dir,
              tag = "states_R848_mult_sepsis_clinSubgrp_pca",
              height = 7,
              width = 7)

# OPLS-DAs
################################################################################

opls_states_LPS <- do_opls_list(imm_states_list_parsed_LPS,
                                slot = "stand",
                                orthoI = 3,
                                permI = 100,
                                plot_dir = paste0(outDir, "opls/states/LPS/plots/"),
                                tag = "LPS")

opls_states_R848 <- do_opls_list(imm_states_list_parsed_R848,
                                 slot = "stand",
                                 orthoI = 3,
                                 permI = 100,
                                 plot_dir = paste0(outDir, "opls/states/R848/plots/"),
                                 tag = "R848")

opls_states_none <- do_opls_list(imm_states_list_parsed_none,
                                 slot = "stand",
                                 orthoI = 3,
                                 permI = 100,
                                 plot_dir = paste0(outDir, "opls/states/none/plots/"),
                                 tag = "none")



opls_states_LPS_sign <- get_opls_sign(opls_states_LPS, vip_thrshld = 1,
                                      save_dir = paste0(outDir,
                                                        "opls/states/LPS/csv/"))
opls_states_R848_sign <- get_opls_sign(opls_states_R848, vip_thrshld = 1,
                                       save_dir = paste0(outDir,
                                                         "opls/states/R848/csv/"))
opls_states_none_sign <- get_opls_sign(opls_states_none, vip_thrshld = 1,
                                       save_dir = paste0(outDir,
                                                         "opls/states/none/csv/"))

# Do OPLS of aid vs sepsis
imm_states_list_parsed_LPS_no_undef <- lapply(imm_states_list_parsed_LPS,
                                              filt_samps,
                                              keep_var = "group",
                                              keep_val = c("aid", "sepsis"))
imm_states_list_parsed_R848_no_undef <- lapply(imm_states_list_parsed_R848,
                                               filt_samps,
                                               keep_var = "group",
                                               keep_val = c("aid", "sepsis"))
imm_states_list_parsed_none_no_undef <- lapply(imm_states_list_parsed_none,
                                               filt_samps,
                                               keep_var = "group",
                                               keep_val = c("aid", "sepsis"))


opls_states_LPS_aid_vs_sepsis <- do_opls_list(imm_states_list_parsed_LPS_no_undef,
                                              slot = "stand",
                                              orthoI = 3,
                                              permI = 100,
                                              plot_dir = paste0(outDir, "opls/states/LPS/plots/"),
                                              tag = "LPS")

get_summary_df(opls_states_LPS_aid_vs_sepsis)

opls_states_R848_aid_vs_sepsis <- do_opls_list(imm_states_list_parsed_R848_no_undef,
                                               slot = "stand",
                                               orthoI = 3,
                                               permI = 100,
                                               plot_dir = paste0(outDir, "opls/states/R848/plots/"),
                                               tag = "R848")

opls_states_none_aid_vs_sepsis <- do_opls_list(imm_states_list_parsed_none_no_undef,
                                               slot = "stand",
                                               orthoI = 3,
                                               permI = 100,
                                               plot_dir = paste0(outDir, "opls/states/none/plots/"),
                                               tag = "none")

opls_summ_stts_LPS_aid_vs_sep <- get_summary_df(opls_states_LPS_aid_vs_sepsis)
opls_summ_stts_R848_aid_vs_sep <- get_summary_df(opls_states_R848_aid_vs_sepsis)
opls_summ_stts_none_aid_vs_sep <- get_summary_df(opls_states_none_aid_vs_sepsis)

opls_summ_stts_LPS_aid_vs_sep$treatment <- "LPS"
opls_summ_stts_R848_aid_vs_sep$treatment <- "R848"
opls_summ_stts_none_aid_vs_sep$treatment <- "none"

opls_summ_stts_aid_vs_sep <- rbind.data.frame(opls_summ_stts_none_aid_vs_sep,
                                              opls_summ_stts_LPS_aid_vs_sep,
                                              opls_summ_stts_R848_aid_vs_sep)
opls_summ_stts_aid_vs_sep <- opls_summ_stts_aid_vs_sep[, c(1:6, 9, 10, 11)]
opls_summ_stts_aid_vs_sep[, c(ncol(opls_summ_stts_aid_vs_sep),
                              1:(ncol(opls_summ_stts_aid_vs_sep) - 1))]

################################################################################
# Best prediction power is in counts l5 none (aid vs sepsis). Let's refine     #
# the model.                                                                   #
################################################################################

opls_RFE <- function(DF,                    # Data
                     y,                     # Response variable.
                     orthoI = 3,            # Number of orthogonal components.
                     worseFeatsProp = NULL) # The proportion of features that
        # are going to be removed on each
        # iteration. 
{
        oplsFit <- opls(DF,
                        y = y,
                        orthoI = orthoI,
                        permI = 2,
                        fig.pdfC = "none")
        r2y <- oplsFit@summaryDF$`R2Y(cum)`
        q2y <- oplsFit@summaryDF$`Q2(cum)`
        
        r2y_vec <- r2y
        q2y_vec <- q2y
        print(dim(DF))
        nFeats <- ncol(DF)
        nFeatsVec <- nFeats
        featsInMod <- paste(colnames(DF), collapse = ", ")
        featsInModVec <- featsInMod
        while(nFeats >= (orthoI + 4)){
                if(is.null(worseFeatsProp)){
                        nWorseFeats <- 1
                }else{
                        nWorseFeats <- ceiling(length(oplsFit@vipVn) * worseFeatsProp)
                }
                # print(nWorseFeats)
                worseFeat <- sort(oplsFit@vipVn, decreasing = F)[1:nWorseFeats]
                print("Features removed:")
                print(paste(names(worseFeat), collapse = ", "))
                DF <- DF[, !colnames(DF) %in% names(worseFeat)]
                print(dim(DF))
                oplsFit <- opls(DF,
                                y = y,
                                orthoI = orthoI,
                                permI = 2,
                                fig.pdfC = "none")
                r2y <- oplsFit@summaryDF$`R2Y(cum)`
                q2y <- oplsFit@summaryDF$`Q2(cum)`
                r2y_vec <- c(r2y_vec, r2y)
                q2y_vec <- c(q2y_vec, q2y)
                nFeats <- ncol(DF)
                nFeatsVec <- c(nFeatsVec, nFeats)
                featsInMod <- paste(colnames(DF), collapse = ", ")
                featsInModVec <- c(featsInModVec, featsInMod)
        }
        statsDF <- data.frame(nFeats = nFeatsVec,
                              r2y = r2y_vec,
                              q2y = q2y_vec,
                              featsInMod = featsInModVec)
        return(statsDF)
}

getBestOplsMod <- function(RFE_DF,
                           protDF,
                           max_n_feats = ncol(protDF),
                           orthoI = 3, respVar = respVar){
        RFE_DF <- RFE_DF[RFE_DF$nFeats <= max_n_feats, ]
        bestProts <- RFE_DF$featsInMod[which.max(RFE_DF$q2y)]
        print(bestProts)
        bestProts <- strsplit(bestProts, split = ", ")[[1]]
        protDF_bestProts <- protDF[, bestProts]
        bestOplsMod <- opls(protDF_bestProts,
                            respVar,
                            orthoI = orthoI,
                            predI = 1,
                            permI = 200)
        return(bestOplsMod)
}


opls_counts_none_aid_vs_sepsis$l5$aid_vs_sepsis@modelDF


counts_l5_none_aid_vs_sepsis_oplsda <- opls(imm_counts_list_parsed_none_no_undef$l5$stand,
                                            y = imm_counts_list_parsed_none_no_undef$l5$sample_info$group,
                                            predI = 1,
                                            orthoI = 2,
                                            permI = 100)

rfe_res <- opls_RFE(imm_counts_list_parsed_none_no_undef$l5$stand,
                    imm_counts_list_parsed_none_no_undef$l5$sample_info$group,
                    orthoI = 2,
                    worseFeatsProp = 0.002)

l5_counts_aid_vs_sepsis_best_mod <- getBestOplsMod(rfe_res,
                                                   imm_counts_list_parsed_none_no_undef$l5$stand,
                                                   respVar = imm_counts_list_parsed_none_no_undef$l5$sample_info$group,
                                                   orthoI = 1,
                                                   max_n_feats = 500)
l5_counts_aid_vs_sepsis_best_mod@modelDF
counts_l5_none_aid_vs_sepsis_oplsda@modelDF

library(caret)
l5_counts_undefined <- filt_samps(imm_counts_list_parsed_none$l5, keep_var = "group", keep_val = "undefined")

l5_counts_4_aid_vs_sepsis_rf <- imm_counts_list_parsed_none_no_undef$l5$log2[, names(l5_counts_aid_vs_sepsis_best_mod@vipVn)]
l5_counts_4_aid_vs_sepsis_rf$group <- imm_counts_list_parsed_none_no_undef$l5$sample_info$group
trControl <- trainControl(method = "cv",
                          number = 10,
                          search = "grid",
                          allowParallel = T)
set.seed(123)
trainFit <- train(group ~ . ,
                  data = l5_counts_4_aid_vs_sepsis_rf,
                  method = "rf",
                  #metric = ,
                  trControl = trControl)

predict(trainFit, l5_counts_undefined$log2)

ctrl <- rfeControl(functions=rfFuncs, method="repeatedcv", number=10)

# With subgroups
l5_counts_4_aid_vs_sepsis_rf <- imm_counts_list_parsed_none_no_undef$l5$log2[, names(l5_counts_aid_vs_sepsis_best_mod@vipVn)]
l5_counts_4_aid_vs_sepsis_rf$group <- imm_counts_list_parsed_none_no_undef$l5$sample_info$subgroup
trControl <- trainControl(method = "cv",
                          number = 10,
                          search = "grid",
                          allowParallel = T)
set.seed(123)
trainFit <- train(group ~ . ,
                  data = l5_counts_4_aid_vs_sepsis_rf,
                  method = "rf",
                  #metric = ,
                  trControl = trControl)

predict(trainFit, l5_counts_undefined$log2)

# RFE with groups
ctrl <- rfeControl(functions=rfFuncs, method="repeatedcv", number=10)

table(imm_counts_list_parsed_none_no_undef$l5$sample_info$group)

set.seed(123)
sizes = c(1:100)
rfe_result <- rfe(x = imm_counts_list_parsed_none_no_undef$l5$log2,
                  y = as.factor(imm_counts_list_parsed_none_no_undef$l5$sample_info$group),
                  sizes=sizes,
                  rfeControl=ctrl)
intersect(rfe_result$optVariables,
          names(l5_counts_aid_vs_sepsis_best_mod@vipVn))

l5_counts_4_aid_vs_sepsis_rf_rfe <- imm_counts_list_parsed_none_no_undef$l5$log2[, rfe_result$optVariables]
l5_counts_4_aid_vs_sepsis_rf_rfe$group <- imm_counts_list_parsed_none_no_undef$l5$sample_info$group

set.seed(123)
trainFit_rfe <- train(group ~ . ,
                      data = l5_counts_4_aid_vs_sepsis_rf_rfe,
                      method = "rf",
                      #metric = ,
                      trControl = trControl)

trainFit_rfe_vImp <- varImp(trainFit_rfe)

trainFit_rfe_vImp_df <- data.frame(feature = rownames(trainFit_rfe_vImp$importance)[order(trainFit_rfe_vImp$importance$Overall, decreasing = T)],
                                   importance = trainFit_rfe_vImp$importance[order(trainFit_rfe_vImp$importance$Overall, decreasing = T), ])

pred_df <- data.frame(sample = rownames(l5_counts_undefined$log2),
                      pred_group = predict(trainFit_rfe, l5_counts_undefined$log2))

# RFE with subgroups

table(imm_counts_list_parsed_none_no_undef$l5$sample_info$subgroup)

set.seed(123)
sizes = c(1:100)
rfe_result <- rfe(x = imm_counts_list_parsed_none_no_undef$l5$log2,
                  y = as.factor(imm_counts_list_parsed_none_no_undef$l5$sample_info$subgroup),
                  sizes=sizes,
                  rfeControl=ctrl)
#intersect(rfe_result$optVariables,
#          names(l5_counts_aid_vs_sepsis_best_mod@vipVn))

l5_counts_4_aid_vs_sepsis_rf_rfe <- imm_counts_list_parsed_none_no_undef$l5$log2[, rfe_result$optVariables]
l5_counts_4_aid_vs_sepsis_rf_rfe$group <- imm_counts_list_parsed_none_no_undef$l5$sample_info$subgroup

set.seed(123)
trainFit_rfe <- train(group ~ . ,
                      data = l5_counts_4_aid_vs_sepsis_rf_rfe,
                      method = "rf",
                      #metric = ,
                      trControl = trControl)

trainFit_rfe_vImp <- varImp(trainFit_rfe)

trainFit_rfe_vImp_df <- data.frame(feature = rownames(trainFit_rfe_vImp$importance)[order(trainFit_rfe_vImp$importance$Overall, decreasing = T)],
                                   importance = trainFit_rfe_vImp$importance[order(trainFit_rfe_vImp$importance$Overall, decreasing = T), ])

pred_df <- data.frame(sample = rownames(l5_counts_undefined$log2),
                      pred_group = predict(trainFit_rfe, l5_counts_undefined$log2))


pred_probs_df <- predict(trainFit_rfe, l5_counts_undefined$log2, type = "prob")
colnames(pred_probs_df) <- paste0("prob_", colnames(pred_probs_df))
pred_df <- cbind.data.frame(pred_df, pred_probs_df)
min(pred_df[pred_df$pred_group == "aid", "prob_aid"])
min(pred_df[pred_df$pred_group == "sepsis_1", "prob_sepsis_1"])
min(pred_df[pred_df$pred_group == "sepsis_2", "prob_sepsis_2"])
min(pred_df[pred_df$pred_group == "sepsis_3", "prob_sepsis_3"])

probs_df_train <- predict(trainFit_rfe, l5_counts_4_aid_vs_sepsis_rf_rfe[, colnames(l5_counts_4_aid_vs_sepsis_rf_rfe) != "group"], type = "prob")
probs_df_train$group <- l5_counts_4_aid_vs_sepsis_rf_rfe$group
min(probs_df_train[probs_df_train$group == "aid", "aid"])
min(probs_df_train[probs_df_train$group == "sepsis_1", "sepsis_1"])
min(probs_df_train[probs_df_train$group == "sepsis_2", "sepsis_2"])
min(probs_df_train[probs_df_train$group == "sepsis_3", "sepsis_3"])

if (!require("entropy", quietly = T)){
        install.packages("entropy")
}
library(entropy)

calc_entropy <- function(votes) {
        probs <- table(votes) / length(votes)  # Convert votes to probabilities
        return(entropy(probs))
}

calc_gini_impurity <- function(votes) {
        probs <- table(votes) / length(votes)  # Convert votes to probabilities
        return(1 - sum(probs^2))
}

pred_trees <- predict(trainFit_rfe$finalModel, l5_counts_undefined$log2, predict.all = TRUE)$individual

pred_df$entropy <- apply(pred_trees, 1, calc_entropy)
pred_df$gini_impurity <- apply(pred_trees, 1, calc_gini_impurity)

pca_undefined <- prcomp(plotUtils::stand(l5_counts_undefined$log2, scale = T, center = T),
                        scale. = F,
                        center = F)
plotUtils::plotPCA(pca_undefined, samp_info = pred_df, col = "gini_impurity",
                   fix_coord = F)
plotUtils::plotPCA(pca_undefined, samp_info = pred_df, col = "pred_group",
                   fix_coord = F)

l5_counts_none_samp_info_pred_undef <- imm_counts_list_parsed_none$l5$sample_info
l5_counts_none_samp_info_pred_undef$subgroup_undef_preds <- l5_counts_none_samp_info_pred_undef$subgroup
l5_counts_none_samp_info_pred_undef$subgroup_undef_preds[l5_counts_none_samp_info_pred_undef$subgroup_undef_preds == "undefined"] <- paste("undefined",
                                                                                                                                           pred_df$pred_group[match(l5_counts_none_samp_info_pred_undef$sample[l5_counts_none_samp_info_pred_undef$subgroup_undef_preds == "undefined"],
                                                                                                                                                                    pred_df$sample)],
                                                                                                                                           sep = "_to_")

plotUtils::plotPCA(PC = imm_counts_pcas_none$l5, samp_info = l5_counts_none_samp_info_pred_undef,
                   col = "subgroup",
                   fix_coord = F)

ggsave(filename = sprintf("%sl5_counts_none_sepsis_subgroups.pdf",
                          outDir),
       height = 4, width = 5)

plotUtils::plotPCA(PC = imm_counts_pcas_none$l5, samp_info = l5_counts_none_samp_info_pred_undef,
                   col = "subgroup_undef_preds",
                   fix_coord = F)

ggsave(filename = sprintf("%sl5_counts_none_sepsis_subgroups_undefined_pred.pdf",
                          outDir),
       height = 4, width = 5)


imm_counts_list_parsed_non








