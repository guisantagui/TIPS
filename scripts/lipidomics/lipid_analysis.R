library(dplyr)

# Create directory if it doesn't exist
createIfNot <- function(pth){
        if(!dir.exists(pth)){
                dir.create(pth, recursive = T)
        }
}

li01_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/data/TIPS-Lipidomics/LI01_Lipidomics_RR/LI01_A1_Data_QS-22-789.xlsx"
li02_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/data/TIPS-Lipidomics/LI02_Lipidomics/LI02_A1_Data_QS-22-789.xlsx"
li03_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/data/TIPS-Lipidomics/LI03_Lipidomics/Lipotype_Report_Wiedemuth_(QS-22-789)(1)/LI03_A1_Data_QS-22-789.xlsx"
li04_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/data/TIPS-Lipidomics/LI04_Lipidomics/Lipotype_Report_Wiedemuth_(QS-24-1124)/A1_Data_QS-24-1124.xlsx"

outDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/results/lipidomics/"

createIfNot(outDir)

clin_dat_file_cat <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/data/TIPS-Medical_data/RedCAP_categorial_values_2024-06-21.csv"
clin_dat_file_num <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/data/TIPS-Medical_data/RedCAP_numerical_values_2024-06-21.csv"

group_data_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/data/TIPS-Data_analysis/TIPS-Results_collection/TIPS_data_collection_2024-07-22/grouping_master_data.RDS"

group_data <- readRDS(group_data_f)

parseLipDat <- function(f){
        #f <- li04_f
        pmol <- data.frame(readxl::read_xlsx(f, sheet = 1))
        molp <- data.frame(readxl::read_xlsx(f, sheet = 2))
        lipDat <- data.frame(readxl::read_xlsx(f, sheet = 3))
        sampDat <- data.frame(readxl::read_xlsx(f, sheet = 4))

        pmol_rn <- lipDat$feature[match(pmol$feature,
                                        lipDat$feature)]

        molp_rn <- lipDat$feature[match(molp$feature,
                                        lipDat$feature)]

        pmol <- pmol[!is.na(pmol_rn), ]
        molp <- molp[!is.na(molp_rn), ]
        rownames(pmol) <- pmol_rn[!is.na(pmol_rn)]
        rownames(molp) <- molp_rn[!is.na(molp_rn)]
        pmol <- pmol[, !colnames(pmol) %in% c("feature", "class")]
        molp <- molp[, !colnames(molp) %in% c("feature", "class")]

        colnames(pmol) <- unlist(as.vector(sampDat[1, match(colnames(pmol),
                                                            colnames(sampDat))]))
        colnames(molp) <- unlist(as.vector(sampDat[1, match(colnames(molp),
                                                            colnames(sampDat))]))

        sampDat <- data.frame(t(sampDat))
        colnames(sampDat) <- sampDat[1, ]
        sampDat <- sampDat[2:nrow(sampDat), ]
        sampDat$shortname <- gsub("X", "", rownames(sampDat))
        rownames(sampDat) <- 1:nrow(sampDat)

        out <- list(pmol = pmol,
                    molp = molp,
                    lipDat = lipDat,
                    sampDat = sampDat)
        return(out)
}

# Load and parse the data
################################################################################
clin_dat_cat <- read.table(clin_dat_file_cat, sep = "\t", header = T)
clin_dat_num <- read.table(clin_dat_file_num, sep = "\t", header = T)

clin_dat_num$Actual.age..years. <- as.numeric(gsub(",",
                                                   ".",
                                                   clin_dat_num$Actual.age..years.))

li01 <- parseLipDat(li01_f)
li02 <- parseLipDat(li02_f)
li03 <- parseLipDat(li03_f)
li04 <- parseLipDat(li04_f)

lipDict <- rbind.data.frame(li01$lipDat,
                            li02$lipDat,
                            li03$lipDat,
                            li04$lipDat)

lipDict <- lipDict[!duplicated(lipDict$feature), ]

# Parse metadata
li01$sampDat$batch <- "li01"
li02$sampDat$batch <- "li02"
li03$sampDat$batch <- "li03"
li04$sampDat$batch <- "li04"

sampColNames <- intersect(colnames(li01$sampDat),
                          colnames(li02$sampDat)) %>%
        intersect(colnames(li03$sampDat)) %>%
        intersect(colnames(li04$sampDat))


sampDatAll <- rbind.data.frame(li01$sampDat[, sampColNames],
                               li02$sampDat[, sampColNames],
                               li03$sampDat[, sampColNames],
                               li04$sampDat[, sampColNames])

li01_toBind <- li01$molp[lipDict$feature, ]
li02_toBind <- li02$molp[lipDict$feature, ]
li03_toBind <- li03$molp[lipDict$feature, ]
li04_toBind <- li04$molp[lipDict$feature, ]

rownames(li01_toBind) <- lipDict$feature
rownames(li02_toBind) <- lipDict$feature
rownames(li03_toBind) <- lipDict$feature
rownames(li04_toBind) <- lipDict$feature

liAll <- cbind.data.frame(li01_toBind,
                          li02_toBind,
                          li03_toBind,
                          li04_toBind)

# Missing value treatment
################################################################################
liAll <- liAll[!apply(liAll, 1, function(x) sum(is.na(x))/length(x) >= .2), ]

liAll <- liAll[, !apply(liAll, 2, function(x) sum(is.na(x))/length(x) >= 0.8)]


library(impute)
liAll <- impute::impute.knn(as.matrix(liAll), rowmax = .2, colmax = .8)$data


sampDatAll$fullname[!sampDatAll$fullname %in% group_data$patient_ID]

liAll <- liAll[, colnames(liAll) %in% group_data$patient_ID[!is.na(group_data$group)]]

liAll_log2 <- t(log2(liAll))

liAll_log2_pca <- prcomp(liAll_log2, scale. = T, center = T)

samp_info <- data.frame(sample = rownames(liAll_log2),
                        group = group_data$group[match(rownames(liAll_log2),
                                                       group_data$patient_ID)],
                        batch = sampDatAll$batch[match(rownames(liAll_log2),
                                                       sampDatAll$fullname)])

library(ggplot2)
library(factoextra)
library(ggrepel)
library(ggpubr)
# Obtaining the top N contributor variables
# to the specified components
getTopContrib <- function(PC, topN = 12, x = "PC1", y = "PC2"){
        if (class(PC) != "prcomp"){
                stop("PC needs to be a prcomp object", call. = F)
        }
        if (topN > nrow(PC$rotation)){
                stop("topN can't be greater than the number of variables in PC object.",
                     call. = F)
        }
        if ((!x %in% colnames(PC$x)) | (!y %in% colnames(PC$x))){
                stop(sprintf("Both %s and %s need to be PCs in PC object",
                             x,
                             y),
                     call. = F)
        }
        contrib <- facto_summarize(PC,
                                   "var",
                                   axes = c(as.numeric(gsub("PC", "", x)),
                                            as.numeric(gsub("PC", "", y))))
        contrib <- contrib[order(contrib$contrib, decreasing = T),
                           c("name", "contrib")]
        topContrib <- as.character(contrib$name[1:topN])
        return(topContrib)
}

getTopContrib(liAll_log2_pca, topN = 10)

dim(liAll_log2_pca$rotation)

plotPCA <- function(PC,
                    x = "PC1",
                    y = "PC2",
                    samp_info = NULL,
                    col = NULL,
                    shape = NULL,
                    labs = F,
                    topNFeats = NULL,
                    biplot = F,
                    fix_coord = T,
                    point_size = 3){
        if (class(PC) != "prcomp"){
                stop("PC needs to be a prcomp object", call. = F)
        }
        dat <- data.frame(obsnames = row.names(PC$x), PC$x)
        if ((!x %in% colnames(dat)[colnames(dat) != "obsnames"]) |
            (!x %in% colnames(dat)[colnames(dat) != "obsnames"])){
                stop(sprintf("%s or %s are not proper principal components.",
                             x,
                             y),
                     call. = F)
        }

        x_sym <- rlang::sym(x)
        y_sym <- rlang::sym(y)
        labs_in <- "obsnames"
        aes_args <- list(x = x_sym, y = y_sym, label = ensym(labs_in))

        dat <- dat[, c("obsnames", x, y)]
        if (!is.null(samp_info)){
                if (class(samp_info) != "data.frame"){
                        stop("sample_info must be a data.frame.",
                             call. = F)
                }
                if (!"sample" %in% colnames(samp_info)){
                        stop("There is no 'sample' column in sample_info dataframe.",
                             call. = F)
                }
                dat <- cbind.data.frame(dat,
                                        samp_info[match(dat$obsnames,
                                                        samp_info$sample),
                                                  colnames(samp_info) != "sample"])
                if (!is.null(col)){
                        if (col %in% colnames(dat)){
                                aes_args$col <- sym(col)
                        }else{
                                stop(sprintf("%s not in sample_info", col),
                                     call. = F)
                        }
                }
                if (!is.null(shape)){
                        if (shape %in% colnames(dat)){
                                aes_args$shape <- sym(shape)
                        }else{
                                stop(sprintf("%s not in sample_info", shape),
                                     call. = F)
                        }
                }
        }
        propVar <- summary(PC)$importance[2, c(x, y)]
        propX <- round(propVar[names(propVar) == x]*100, digits = 2)
        propY <- round(propVar[names(propVar) == y]*100, digits = 2)

        pcaPlt <- ggplot(dat, do.call(aes, aes_args)) +
                geom_point(size = point_size) +
                xlab(sprintf("%s (%s %%)", x, propX)) +
                ylab(sprintf("%s (%s %%)", y, propY)) +
                theme(title = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA,
                                                  linewidth = 1),
                      panel.grid.major = element_line(colour = "#d4d4d4"),
                      legend.position = "right")
        if (fix_coord){
                pcaPlt <- pcaPlt +
                        coord_fixed()
        }
        if (labs){
                pcaPlt <- pcaPlt +
                        geom_text_repel()
        }
        if (biplot){
                datapc <- data.frame(varnames=rownames(PC$rotation),
                                     PC$rotation)
                mult <- min(
                        (max(dat[,y]) - min(dat[,y])/(max(datapc[,y])-min(datapc[,y]))),
                        (max(dat[,x]) - min(dat[,x])/(max(datapc[,x])-min(datapc[,x])))
                )
                datapc <- transform(datapc,
                                    v1 = .7 * mult * (get(x)),
                                    v2 = .7 * mult * (get(y))
                )
                datapc$x0 <- rep(0, nrow(datapc))
                datapc$y0 <- rep(0, nrow(datapc))
                if(!is.null(topNFeats)){
                        varPlotFilt <- getTopContrib(PC,
                                                     topN = topNFeats,
                                                     x = x,
                                                     y = y)
                        datapc <- datapc[datapc$varnames %in% varPlotFilt, ]
                }
                pcaPlt <- pcaPlt +
                        geom_text_repel(data=datapc,
                                        aes(x=v1, y=v2,
                                            label=varnames),
                                        color = "black",
                                        size = 3,
                                        max.overlaps = 100,
                                        inherit.aes = F) +
                        geom_segment(data = datapc, aes(x=x0,
                                                        y=y0,
                                                        xend=v1,
                                                        yend=v2),
                                     arrow = arrow(length=unit(0.2,"cm"),
                                                   type = "closed",
                                                   angle = 20),
                                     alpha=0.75,
                                     color="black",
                                     linewidth = 0.5,
                                     inherit.aes = F)
        }
        return(pcaPlt)
}

doPCAMultiPlot <- function(PC,
                           nComps,
                           samp_info = NULL,
                           col = NULL,
                           shape = NULL,
                           labs = F,
                           topNFeats = NULL,
                           biplot = F,
                           point_size = 3){
        if (class(PC)[1] != "prcomp"){
                stop("PC needs to be a prcomp object", call. = F)
        }
        plotList <- list()
        for(j in 2:(nComps)){
                #j <- 2
                for(i in 1:(nComps - 1)){
                        #i <- 1
                        if(j > i){
                                scPlot <- plotPCA(PC,
                                                  x = sprintf("PC%s", i),
                                                  y = sprintf("PC%s", j),
                                                  samp_info = samp_info,
                                                  col = col,
                                                  shape = shape,
                                                  labs = labs,
                                                  topNFeats = topNFeats,
                                                  biplot = biplot,
                                                  fix_coord = F,
                                                  point_size = point_size)
                                if(j < nComps){
                                        scPlot <- scPlot +
                                                theme(axis.title.x = element_blank(),
                                                      axis.text.x = element_blank())
                                }
                                if(i > 1){
                                        scPlot <- scPlot +
                                                theme(axis.title.y = element_blank(),
                                                      axis.text.y = element_blank())
                                }
                        }else{
                                scPlot <- NA
                        }
                        plotList[[sprintf("PC%s_PC%s", i, j)]] <- scPlot
                }
        }
        multPlot <- ggarrange(plotlist = plotList,
                              common.legend = T,
                              ncol = nComps - 1,
                              nrow = nComps - 1,
                              widths = c(1, rep(.8, nComps-2)),
                              heights = c(rep(.8, nComps-2), 1))
        return(multPlot)
}

plt <- plotPCA(liAll_log2_pca,
        x = "PC1",
        y = "PC2",
        samp_info = samp_info,
        col = "group",
        topNFeats = 10,
        biplot = T,
        fix_coord = F)

pca_mult_plot_batch <- doPCAMultiPlot(liAll_log2_pca, nComps = 5,
                                      samp_info = samp_info,
                                      col = "batch",
                                      topNFeats = 5,
                                      biplot = T,
                                      point_size = 1)

ggsave(filename = sprintf("%spca_mult_plot_batch.pdf", outDir),
       plot = pca_mult_plot_batch)

pca_mult_plot_pheno <- doPCAMultiPlot(liAll_log2_pca, nComps = 5,
                                      samp_info = samp_info,
                                      col = "group",
                                      topNFeats = 5,
                                      biplot = T,
                                      point_size = 1)

ggsave(filename = sprintf("%spca_mult_plot_pheno.pdf", outDir),
       plot = pca_mult_plot_pheno, width = 10, height = 7.5)

library(ropls)


y_aid_vs_rest <- samp_info$group[match(rownames(liAll_log2), samp_info$sample)]
y_aid_vs_rest[y_aid_vs_rest != "aid"] <- "rest"
pdf(file = sprintf("%soplsda_aid_vs_rest.pdf", outDir))
oplsda_aid_vs_rest <- opls(liAll_log2,
                           y = y_aid_vs_rest,
                           predI = 1,
                           orthoI = 4,
                           permI = 200)
dev.off()

y_sepsis_vs_rest <- samp_info$group[match(rownames(liAll_log2), samp_info$sample)]
y_sepsis_vs_rest[y_sepsis_vs_rest != "sepsis"] <- "rest"
pdf(file = sprintf("%soplsda_sepsis_vs_rest.pdf", outDir))
oplsda_sepsis_vs_rest <- opls(liAll_log2,
                              y = y_sepsis_vs_rest,
                              predI = 1,
                              orthoI = 4,
                              permI = 200)
dev.off()

y_undef_vs_rest <- samp_info$group[match(rownames(liAll_log2), samp_info$sample)]
y_undef_vs_rest[y_undef_vs_rest != "undefined"] <- "rest"
pdf(file = sprintf("%soplsda_undef_vs_rest.pdf", outDir))
oplsda_undef_vs_rest <- opls(liAll_log2,
                             y = y_undef_vs_rest,
                             predI = 1,
                             orthoI = 4,
                             permI = 200)
dev.off()

oplsda_aid_vs_rest_sign <- oplsda_aid_vs_rest@vipVn[oplsda_aid_vs_rest@vipVn >= 1]
oplsda_sepsis_vs_rest_sign <- oplsda_sepsis_vs_rest@vipVn[oplsda_sepsis_vs_rest@vipVn >= 1]

# Given a oplsda object and 
get_opls_lip_cat_bplot <- function(opls_obj, thrshld = 1, var_info){
        opls_sign <- opls_obj@vipVn[opls_obj@vipVn >= thrshld]
        cats <- var_info$class[match(names(opls_sign), var_info$feature)]
        cat_tab <- table(cats)
        dat_df <- data.frame(category = names(cat_tab),
                             frequency = as.vector(cat_tab)/sum(as.vector(cat_tab)))
        plt <- ggplot(dat_df, mapping = aes(x = category,
                                            y = frequency)) +
                geom_bar(stat = "identity") +
                theme(axis.title.x = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA,
                                                  linewidth = 1),
                      panel.grid.major = element_blank(),
                      legend.position = "right")
        return(plt)
}

aid_vs_rest_sign_lipFreq <- get_opls_lip_cat_bplot(oplsda_aid_vs_rest,
                                                   thrshld = 1,
                                                   var_info = lipDict)

ggsave(filename = sprintf("%said_vs_rest_sign_lipFreq.pdf",
                          outDir),
       plot = aid_vs_rest_sign_lipFreq,
       width = 5, height = 3)

global_lipFreq <- get_opls_lip_cat_bplot(oplsda_aid_vs_rest,
                                         thrshld = 0,
                                         var_info = lipDict)

ggsave(filename = sprintf("%sglobal_lipFreq.pdf",
                          outDir),
       plot = global_lipFreq,
       width = 5, height = 3)


get_opls_lip_cat_bplot(oplsda_sepsis_vs_rest, thrshld = 1, var_info = lipDict)


do_cat_ORA <- function(opls_obj, thrshld = 1, var_info, adj_method = "BH"){
        #opls_obj <- oplsda_aid_vs_rest
        #thrshld <- 1
        #var_info <- lipDict
        opls_sign <- names(opls_obj@vipVn[opls_obj@vipVn >= thrshld])
        all_vars <- names(opls_obj@vipVn)
        all_vars_df <- data.frame(var = all_vars,
                                  class = var_info$class[match(all_vars,
                                                               var_info$feature)])
        uniq_cats <- unique(all_vars_df$class)
        p_vals <- c()
        for(i in seq_along(uniq_cats)){
                #i <- 1
                ct <- uniq_cats[i]
                sign_isCat <- sum(all_vars_df$class[match(opls_sign, all_vars_df$var)] == ct)
                sign_isNotCat <- sum(all_vars_df$class[match(opls_sign, all_vars_df$var)] != ct)
                notSign_isCat <- sum(all_vars_df$class == ct) - sign_isCat
                notSign_isNotCat <- sum(all_vars_df$class != ct) - sign_isNotCat
                cont_mat <- matrix(c(sign_isCat, sign_isNotCat,
                                     notSign_isCat, notSign_isNotCat),
                                   nrow = 2, byrow = T,
                                   dimnames = list(c("significant",
                                                    "not_significant"),
                                                   c("in_class",
                                                     "not_in_class")))
                p_val <- fisher.test(cont_mat, alternative = "greater")$p.value
                p_vals <- c(p_vals, p_val)
        }
        res_df <- data.frame(class = uniq_cats,
                             p_val = p_vals,
                             p_adj = p.adjust(p_vals, method = adj_method))
        return(res_df)
}

aid_vs_rest_cat_ora <- do_cat_ORA(oplsda_aid_vs_rest, thrshld = 1, var_info = lipDict)

enzyme_file <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/data/enzymes.tsv.gz"
lipids_file <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/data/lipids.tsv.gz"
lip2up_file <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/data/lipids2uniprot.tsv.gz"
go_file <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/data/go.tsv"

lipids <- data.frame(data.table::fread(lipids_file, sep = "\t", fill = T))
lip2up <- data.frame(data.table::fread(lip2up_file, sep = "\t", fill = T))
enzymes <- data.frame(data.table::fread(enzyme_file, sep = "\t", fill = T))
enzymes[enzymes$Taxon.scientific.name == "Homo sapiens", ]
go_lip <- data.frame(data.table::fread(go_file, sep = "\t", fill = T))
lipDict$SwissLipids.ID[match(colnames(liAll_log2), lipDict$feature)] %in% enzymes$SwissLipids.ID

View(lipids[match(lipDict$SwissLipids.ID[match(colnames(liAll_log2), lipDict$feature)],
                  lipids$Lipid.ID), ])
lipDict$SwissLipids.ID[match(colnames(liAll_log2), lipDict$feature)] %in% lipids$Lipid.ID
lipDict$SwissLipids.ID[match(colnames(liAll_log2), lipDict$feature)] %in% lip2up$metabolite.id
lipDict$SwissLipids.ID[match(colnames(liAll_log2), lipDict$feature)] %in% go_lip$Lipid.ID
go_lip[match(lipDict$SwissLipids.ID[match(names(oplsda_aid_vs_rest_sign), lipDict$feature)][lipDict$SwissLipids.ID[match(names(oplsda_aid_vs_rest_sign), lipDict$feature)] %in% go_lip$Lipid.ID],
             go_lip$Lipid.ID), ]
View(lip2up[match(lipDict$SwissLipids.ID[match(colnames(liAll_log2), lipDict$feature)],
             lip2up$metabolite.id), ])


sum(lipDict$SwissLipids.ID[match(names(oplsda_aid_vs_rest_sign), lipDict$feature)] %in% lip2up$metabolite.id)/length(oplsda_aid_vs_rest_sign)

aid_vs_rest_sign_prots <- lip2up$UniprotKB.IDs[match(lipDict$SwissLipids.ID[match(names(oplsda_aid_vs_rest_sign), lipDict$feature)][lipDict$SwissLipids.ID[match(names(oplsda_aid_vs_rest_sign), lipDict$feature)] %in% lip2up$metabolite.id],
                                                     lip2up$metabolite.id)]

unique(unlist(sapply(aid_vs_rest_sign_prots, function(x) strsplit(x, split = " | ", fixed = T)[[1]])))

library(clusterProfiler)
sort(lipDict$SwissLipids.ID[match(colnames(liAll_log2), lipDict$feature)])

# Generate a txt file with all lipid IDs for use as background in LION.
all_lips <- data.frame(lipID = lipDict$Shorthand.Notation[(match(names(oplsda_aid_vs_rest@vipVn),
                                                                 lipDict$feature))])
write.table(all_lips, file = sprintf("%sbg_lips.txt", outDir),
            quote = F, col.names = F, row.names = F)

aid_vs_rest_signLips <- data.frame(lipID = lipDict$Shorthand.Notation[match(names(oplsda_aid_vs_rest_sign),
                                                                            lipDict$feature)])
write.table(aid_vs_rest_signLips, file = sprintf("%said_vs_rest_signLips.txt", outDir),
            quote = F, col.names = F, row.names = F)

aid_vs_rest_VIPs <- data.frame(lipID = lipDict$Shorthand.Notation[match(names(oplsda_aid_vs_rest@vipVn),
                                                                        lipDict$feature)],
                               VIP = oplsda_aid_vs_rest@vipVn)
write.table(aid_vs_rest_VIPs, file = sprintf("%said_vs_rest_VIPs.txt", outDir),
            quote = F, col.names = F, row.names = F, sep = "\t")

# Load LION enrichment and plot it:
lion_enr_aid_vs_rest_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/results/lipidomics/LION_aid_vs_rest_sign.csv"

lion_enr_aid_vs_rest <- read.csv(lion_enr_aid_vs_rest_f)
lion_enr_aid_vs_rest_sign <- lion_enr_aid_vs_rest[lion_enr_aid_vs_rest$FDR.q.value <= 0.05, ]

get_lion_dotplot <- function(lion, thrshld = 0.05){
        #lion <- lion_enr_aid_vs_rest
        lion <- lion[lion$FDR.q.value <= thrshld, ]
        lion$ratio <- lion$Significant/lion$Annotated
        plt <- ggplot(lion, mapping = aes(x = Discription,
                                          y = ratio,
                                          color = FDR.q.value,
                                          size = Significant)) +
                geom_point() + 
                scale_size(range = c(2, 10)) +
                labs(x = "Description", y = "enrichment ratio",
                     color = "FDR",
                     size = "# of Significant DALs") +
                scale_color_gradient(low = "red", high = "blue") +
                guides(fill = guide_colorbar(reverse = TRUE)) +
                coord_flip() +
                scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +
                theme(axis.text.y = element_text(size=15),
                      axis.text.x = element_text(size=15),
                      axis.title.x = element_text(size=18),
                      axis.title.y = element_blank(),
                      title = element_text(size=20),
                      panel.background = element_blank(),
                      panel.grid.major = element_line(colour = "gray"),
                      panel.grid.minor = element_blank(),
                      legend.text = element_text(size=12),
                      legend.title = element_text(size=13),
                      axis.line = element_line(colour = "black"),
                      axis.line.y = element_line(colour = "black"),
                      panel.border = element_rect(colour = "black",
                                                  fill=NA,
                                                  linewidth = 1))
        return(plt)
}
lion_enr_aid_vs_rest_dtplt <- get_lion_dotplot(lion_enr_aid_vs_rest,
                                               thrshld = 0.05)
ggsave(filename = sprintf("%slion_enr_aid_vs_rest_dtplt.pdf",
                          outDir),
       lion_enr_aid_vs_rest_dtplt, height = 4, width = 10)

lion_db_file <- "/Users/guillem.santamaria/Documents/postdoc/comput/TIPS/data/LION.csv.gz"
lion_db <- data.frame(data.table::fread(lion_db_file))

View(lion_db)
lion_db[match(lion_enr_aid_vs_rest_sign$Term.ID, lion_db$http...www.w3.org.2004.02.skos.core.notation), ]


oplsda_aid_vs_rest_loads <- data.frame(feat = rownames(oplsda_aid_vs_rest@loadingMN),
                                       loading = as.vector(oplsda_aid_vs_rest@loadingMN))
oplsda_aid_vs_rest_loads <- oplsda_aid_vs_rest_loads[order(oplsda_aid_vs_rest_loads$loading), ]
oplsda_aid_vs_rest_loads$feat <- factor(oplsda_aid_vs_rest_loads$feat,
                                        levels = oplsda_aid_vs_rest_loads$feat)

loads_aid_vs_rest_TAG <- loads_aid_vs_rest_TAG[grep("TAG", loads_aid_vs_rest_TAG$feat), ]
loads_aid_vs_rest_TAG <- loads_aid_vs_rest_TAG[order(loads_aid_vs_rest_TAG$loading), ]
loads_aid_vs_rest_TAG$feat <- factor(loads_aid_vs_rest_TAG$feat,
                                     levels = loads_aid_vs_rest_TAG$feat)
loads_aid_vs_rest_TAG_bPlt <- ggplot(data = loads_aid_vs_rest_TAG,
                                     mapping = aes(x = loading, y = feat)) +
        geom_bar(stat = "identity") +
        xlab("Predictive component loading (rest vs aid)") +
        theme(axis.text.y = element_text(size=15),
              axis.text.x = element_text(size=15),
              axis.title.x = element_text(size=18),
              axis.title.y = element_blank(),
              title = element_text(size=20),
              panel.background = element_blank(),
              panel.grid.major = element_line(colour = "gray"),
              panel.grid.minor = element_blank(),
              legend.text = element_text(size=12),
              legend.title = element_text(size=13),
              axis.line = element_line(colour = "black"),
              axis.line.y = element_line(colour = "black"),
              panel.border = element_rect(colour = "black",
                                          fill=NA,
                                          linewidth = 1))
ggsave(filename = sprintf("%sloads_aid_vs_rest_TAG_bPlt.pdf", outDir), plot = loads_aid_vs_rest_TAG_bPlt)

oplsda_aid_vs_rest_loads$loading
martijnmolenaar/topOnto.LION2.db/
devtools::install_github("martijnmolenaar/topOnto.LION2.db", upgrade = "never")
devtools::install_github("martijnmolenaar/topOnto.LION2.db/topOnto.LION.db",
                         upgrade = "never")

if(!require("topOnto.LION.db", quietly = T)){
        devtools::install_github("martijnmolenaar/topOnto.LION2.db", upgrade = "never")
}
if(!require("topOnto", quietly = T)){
        devtools::install_github("hxin/topOnto", upgrade = "never")
}
library(topOnto)
library(topOnto.LION.db)
data(ONTdata)
?topOnto.LION.db::ONT

