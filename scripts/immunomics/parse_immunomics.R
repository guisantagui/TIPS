################################################################################
# TIPS: Parse immunomics files and save them as RDS objects.                   #
################################################################################

if (!require("plotUtils", quietly = T)){
        devtools::install_github("guisantagui/plotUtils", upgrade = "never")
}
library(plotUtils)

# Direcotry stuff
################################################################################
immu_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/data/Immunomics_master_data_2025-01-10.RData"
group_data_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/data/TIPS-Data_analysis/TIPS-Results_collection/TIPS_data_collection_2025-01-11/grouping_master_data_2025-01-10.RDS"
subgroup_data_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/data/TIPS patient classification_sb_cs27022025.xlsx"
outDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/results/immunomics/parsed/"

create_dir_if_not(outDir)

# Functions
################################################################################

parse_immun <- function(immun, patient_info = NULL, type = "counts", max_zero_prop_feats = 0.8){
        #immun <- immunomics_counts_level5
        #patient_info <- group_data
        immun_comm <- data.frame(immun[, 1:grep("dalama", colnames(immun))])
        if (type == "counts"){
                immun_dat <- data.frame(immun[, (grep("dalama",
                                                      colnames(immun)) + 2):ncol(immun)])
                # Remove columns that have high proportion of zero counts
                zero_prop_feats <- apply(immun_dat,
                                         2,
                                         function(x) sum(x == 0) / length(x)) 
                excluded_feats <- names(zero_prop_feats)[zero_prop_feats >= max_zero_prop_feats]
                if(length(excluded_feats) > 0){
                        warning(sprintf("Removing %s columns from data due to having more than %s%% of zeros.",
                                        length(excluded_feats),
                                        max_zero_prop_feats * 100),
                                call. = F)
                }
                immun_dat <- immun_dat[, !colnames(immun_dat) %in% excluded_feats]
                immun_dat_norm <- (immun_dat + 1)/(rowSums(immun_dat) + ncol(immun_dat))
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
                immun_comm$clin_regrouping <- patient_info$clin_regrouping[match(immun_comm$patient_ID,
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

# Get an object of combined treatments, with differential response ratio (DRR),
# (LPS - none)/(R848 - none)
get_drr_obj <- function(pars_obj, slot = "norm"){
        drr_obj <- list()
        for(i in seq_along(pars_obj)){
                lvl <- names(pars_obj)[i]
                none <- filt_samps(pars_obj[[i]],
                                   keep_var = "treatment",
                                   keep_val = "none")
                lps <- filt_samps(pars_obj[[i]],
                                  keep_var = "treatment",
                                  keep_val = "LPS")
                R848 <- filt_samps(pars_obj[[i]],
                                   keep_var = "treatment",
                                   keep_val = "R848")
                comm_samps <- intersect(intersect(none$sample_info$patient_ID,
                                                  lps$sample_info$patient_ID),
                                        R848$sample_info$patient_ID)
                samp_info <- none$sample_info[none$sample_info$patient_ID %in% comm_samps, ]
                samp_info$sample <- gsub("\\_.*", "", samp_info$sample)
                samp_info$treatment <- "drr"
                none <- none[[slot]]
                lps <- lps[[slot]]
                R848 <- R848[[slot]]
                rownames(none) <- gsub("\\_.*", "", rownames(none))
                rownames(lps) <- gsub("\\_.*", "", rownames(lps))
                rownames(R848) <- gsub("\\_.*", "", rownames(R848))
                comm_cols <- intersect(intersect(colnames(none),
                                                 colnames(lps)),
                                       colnames(R848))
                none <- none[match(samp_info$sample, rownames(none)),
                             match(comm_cols, colnames(none))]
                lps <- lps[match(samp_info$sample, rownames(lps)),
                           match(comm_cols, colnames(lps))]
                R848 <- R848[match(samp_info$sample, rownames(R848)),
                             match(comm_cols, colnames(R848))]
                drr <- ((lps - none) + .0000001)/((R848 - none) + .0000001)
                #drr <- (lps)/(R848)
                obj <- list(sample_info = samp_info)
                obj[[slot]] <- drr
                drr_obj[[lvl]] <- obj
        }
        return(drr_obj)
}

# Get difference between two treatment objects. Intended to get the differences
# of either LPS and none, or R848 and none. objF[[slot]] - objI[[slot]]
get_treatcomp_obj <- function(pars_obj, slot = "log2", objF = "LPS", objI = "none"){
        #pars_obj <- imm_counts_list_parsed
        #slot = "log2"
        #objF = "LPS"
        #objI = "none"
        comp_obj <- list()
        for(i in seq_along(pars_obj)){
                lvl <- names(pars_obj)[i]
                objI_lvl <- filt_samps(pars_obj[[i]],
                                       keep_var = "treatment",
                                       keep_val = objI)
                objF_lvl <- filt_samps(pars_obj[[i]],
                                       keep_var = "treatment",
                                       keep_val = objF)
                if (!slot %in% names(objI_lvl)){
                        stop(sprintf("%s is not a slot of pars_obj", slot),
                             call. = F)
                }
                comm_samps <- intersect(objI_lvl$sample_info$patient_ID,
                                        objF_lvl$sample_info$patient_ID)
                samp_info <- objI_lvl$sample_info[objI_lvl$sample_info$patient_ID %in% comm_samps, ]
                samp_info$sample <- gsub("\\_.*", "", samp_info$sample)
                samp_info$treatment <- sprintf("%sminus%s", objF, objI)
                objI_lvl <- objI_lvl[[slot]]
                objF_lvl <- objF_lvl[[slot]]
                rownames(objI_lvl) <- gsub("\\_.*", "", rownames(objI_lvl))
                rownames(objF_lvl) <- gsub("\\_.*", "", rownames(objF_lvl))
                comm_cols <- intersect(colnames(objI_lvl),
                                       colnames(objI_lvl))
                objI_lvl <- objI_lvl[match(samp_info$sample, rownames(objI_lvl)),
                                     match(comm_cols, colnames(objI_lvl))]
                objF_lvl <- objF_lvl[match(samp_info$sample, rownames(objF_lvl)),
                                     match(comm_cols, colnames(objF_lvl))]
                comp <- objF_lvl - objI_lvl
                obj <- list(sample_info = samp_info)
                obj[[slot]] <- comp
                comp_obj[[lvl]] <- obj
        }
        return(comp_obj)
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
        slots <- slots[!grepl("OLs", slots)]
        for (s in slots){
                obj[[s]] <- obj[[s]][rownames(obj[[s]]) %in% obj$sample_info$sample, ]
        }
        return(obj)
}

# Imputes values with KNN over a list, in a given slot
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
# Load data
################################################################################

load(immu_f)
group_data <- readRDS(group_data_f)
subgroup_data <- as.data.frame(readxl::read_xlsx(subgroup_data_f))
subgroup_data <- subgroup_data[!is.na(subgroup_data$`Study id`), ]
group_data$clin_subgroup <- subgroup_data$subgroup[match(group_data$patient_ID, subgroup_data$`Study id`)]
group_data$clin_regrouping <- subgroup_data$regrouping[match(group_data$patient_ID, subgroup_data$`Study id`)]
group_data$clin_subgroup <- gsub(" ", "_", group_data$clin_subgroup)
group_data$clin_subgroup <- gsub("group1", "group_1", group_data$clin_subgroup)
group_data$clin_subgroup <- gsub("Group_2", "group_2", group_data$clin_subgroup)

table(group_data$clin_subgroup[group_data$group == "undefined"])
group_data[group_data$group == "undefined", ]
table(group_data$clin_regrouping[group_data$group == "undefined"])

not_in_groupDat <- unique(immunomics_counts_level1$patient_ID[!immunomics_counts_level1$patient_ID %in% group_data$patient_ID])

# Parse counts files
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

# Plot distributions
lapply(imm_counts_list_parsed, function(x) plot(density(as.matrix(x$stand))))

# Get DRR ((lps - none)/(R848 - none))
imm_counts_list_parsed_drr <- get_drr_obj(imm_counts_list_parsed, slot = "log2")
#imm_counts_list_parsed_drr <- lapply(imm_counts_list_parsed_drr,
#                                     do_log2,
#                                     slot = "norm")
imm_counts_list_parsed_drr <- lapply(imm_counts_list_parsed_drr,
                                     do_stand,
                                     slot = "log2")

# Plot distributions
lapply(imm_counts_list_parsed_drr, function(x) plot(density(as.matrix(x$stand))))

# Get object of treatments minus none
imm_counts_list_parsed_LPSminusNone <- get_treatcomp_obj(imm_counts_list_parsed,
                                                         slot = "log2",
                                                         objF = "LPS",
                                                         objI = "none")

imm_counts_list_parsed_LPSminusNone <- lapply(imm_counts_list_parsed_LPSminusNone,
                                              do_stand,
                                              slot = "log2")

imm_counts_list_parsed_R848minusNone <- get_treatcomp_obj(imm_counts_list_parsed,
                                                          slot = "log2",
                                                          objF = "R848",
                                                          objI = "none")

imm_counts_list_parsed_R848minusNone <- lapply(imm_counts_list_parsed_R848minusNone,
                                               do_stand,
                                               slot = "log2")

# Plot distributions
lapply(imm_counts_list_parsed_LPSminusNone, function(x) plot(density(as.matrix(x$stand))))
lapply(imm_counts_list_parsed_R848minusNone, function(x) plot(density(as.matrix(x$stand))))

# Get objects by treatment
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

# Plot distributions
lapply(imm_counts_list_parsed_none, function(x) plot(density(as.matrix(x$stand))))
lapply(imm_counts_list_parsed_LPS, function(x) plot(density(as.matrix(x$stand))))
lapply(imm_counts_list_parsed_R848, function(x) plot(density(as.matrix(x$stand))))


# Save objects
saveRDS(imm_counts_list_parsed,
        file = sprintf("%simm_counts_parsed.rds",
                       outDir))

saveRDS(imm_counts_list_parsed_none,
        file = sprintf("%simm_counts_parsed_none.rds",
                       outDir))

saveRDS(imm_counts_list_parsed_LPS,
        file = sprintf("%simm_counts_parsed_LPS.rds",
                       outDir))

saveRDS(imm_counts_list_parsed_R848,
        file = sprintf("%simm_counts_parsed_R848.rds",
                       outDir))

saveRDS(imm_counts_list_parsed_drr,
        file = sprintf("%simm_counts_parsed_drr.rds",
                       outDir))

saveRDS(imm_counts_list_parsed_LPSminusNone,
        file = sprintf("%simm_counts_parsed_LPSminusNone.rds",
                       outDir))

saveRDS(imm_counts_list_parsed_R848minusNone,
        file = sprintf("%simm_counts_parsed_R848minusNone.rds",
                       outDir))
# Parse states files
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

lapply(imm_counts_list_parsed, function(x) plot(density(as.matrix(x$stand))))

# Get DRR ((lps - none)/(R848 - none))
imm_states_list_parsed_drr <- get_drr_obj(imm_states_list_parsed, slot = "log2")
#imm_states_list_parsed_drr <- lapply(imm_states_list_parsed_drr,
#                                     do_log2,
#                                     slot = "log2")
imm_states_list_parsed_drr <- lapply(imm_states_list_parsed_drr,
                                     do_stand,
                                     slot = "log2")

# Plot distributions
lapply(imm_states_list_parsed_drr, function(x) plot(density(as.matrix(x$stand))))

# Get objects by treatment
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

lapply(imm_states_list_parsed_none, function(x) plot(density(as.matrix(x$stand))))
lapply(imm_states_list_parsed_LPS, function(x) plot(density(as.matrix(x$stand))))
lapply(imm_states_list_parsed_R848, function(x) plot(density(as.matrix(x$stand))))

# Save objects
saveRDS(imm_states_list_parsed,
        file = sprintf("%simm_states_parsed.rds",
                       outDir))

saveRDS(imm_states_list_parsed_none,
        file = sprintf("%simm_states_parsed_none.rds",
                       outDir))

saveRDS(imm_states_list_parsed_LPS,
        file = sprintf("%simm_states_parsed_LPS.rds",
                       outDir))

saveRDS(imm_states_list_parsed_R848,
        file = sprintf("%simm_states_parsed_R848.rds",
                       outDir))

saveRDS(imm_states_list_parsed_drr,
        file = sprintf("%simm_states_parsed_drr.rds",
                       outDir))
