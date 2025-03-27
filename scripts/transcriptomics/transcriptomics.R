if(!require("nanostringr", quietly = T)){
        BiocManager::install("nanostringr", update = F)
}
library(nanostringr)
library(plotUtils)

data_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/data/TIPS-Data_analysis/TIPS-Results_collection/TIPS_data_collection_2024-07-22/Transcriptomics_master_data.RDS"
group_data_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/data/TIPS-Data_analysis/TIPS-Results_collection/TIPS_data_collection_2025-01-11/grouping_master_data_2025-01-10.RDS"
rcc <- read_rcc("/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/data/TIPS-Tanscriptomics/Experiment_attachments_12_2023___Transcriptomics_Batch_1/Files/NanoString_raw data/VOL16_20231128_210362720724_RCC/")

dim(rcc$raw)

d <- readRDS(data_f)
group_data <- readRDS(group_data_f)
colnames(group_data) <- gsub("patient_ID", "sample", colnames(group_data))
d_filt <- d[, !grepl("NEG|POS", colnames(d))]

View(d_filt)
rownames(d_filt) <- d_filt$patient_ID

run_info <- d_filt[, 1:7]
colnames(run_info) <- gsub("patient_ID", "sample", colnames(run_info))
d_filt <- d_filt[, 8:ncol(d_filt)]

d_filt <- as.data.frame(apply(d_filt, 2, function(x) as.numeric(gsub(",", ".", x))))
rownames(d_filt) <- rownames(run_info)



d_filt_log <- log2(d_filt)
d_filt_log_stand <- stand(d_filt_log)
d_filt_log_stand_combat <- t(sva::ComBat(t(d_filt_log_stand),
                                         batch = run_info$batch_experiment))

plot(density(as.matrix(d_filt_log_stand)))

d_filt_log_stand_pca <- prcomp(d_filt_log_stand, scale. = F, center = F)
d_filt_log_stand_combat_pca <- prcomp(d_filt_log_stand_combat, scale. = F, center = F)

rownames(d_filt_log_stand)[!rownames(d_filt_log_stand) %in% group_data$sample]
plotPCA(d_filt_log_stand_pca, samp_info = group_data, col = "group")
plotPCA(d_filt_log_stand_pca, samp_info = run_info, col = "batch_experiment")

plotPCA(d_filt_log_stand_combat_pca, samp_info = group_data, col = "group")

doPCAMultiPlot(d_filt_log_stand_combat_pca,
               samp_info = group_data,
               col = "group",
               nComps = 5)
plotPCA(d_filt_log_stand_combat_pca,
        samp_info = run_info,
        col = "batch_experiment")

plotPCA(d_filt_log_stand_pca, samp_info = run_info, col = "batch_experiment")
