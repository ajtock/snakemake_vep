#!/usr/local/bin/Rscript

# Author: Andy Tock
# Date: 13/12/22

# Usage:
# ./cluster_analysis.R results/ vaf_DF_dp_patho_merge_vaf_DF_dp.tsv
# ./cluster_analysis.R results/ vaf_DF_dp_patho.tsv
# ./cluster_analysis.R results/ vaf_DF_dp_CRC_Control_patho_merge_vaf_DF_dp_CRC_Control.tsv
# ./cluster_analysis.R results/ vaf_DF_dp_CRC_Control_patho.tsv

#indir="results/"
#vaf_tsv="vaf_DF_dp_patho_merge_vaf_DF_dp.tsv"

args = commandArgs(trailingOnly=T)
indir = args[1]
vaf_tsv = args[2]

options(stringsAsFactors=F)
library(tidyr)
library(dplyr)
library(ComplexHeatmap)
library(magick)
library(circlize)
library(RColorBrewer)
library(viridis)
library(data.table)
library(fpc)
library(cluster)

outdir = paste0(indir, "cluster_analysis/")
plotdir = paste0(outdir, "plots/")
system(paste0("[ -d ", outdir, " ] || mkdir -p ", outdir))
system(paste0("[ -d ", plotdir, " ] || mkdir -p ", plotdir))

# Load VAF table as data.frame
full_DF = read.table(paste0(indir, vaf_tsv), 
                     header=F,
                     sep="\t",
                     na.strings="NaN")
#full_DF = fread(paste0(indir, vaf_tsv),
#                header=F,
#                sep="\t",
#                na.strings="NaN",
#                fill=F,
#                quote="\"",
#                check.names=F,
#                data.table=F)
# Get VAF matrix
vaf_DF = full_DF[6:nrow(full_DF), 7:ncol(full_DF)]
vaf_DF[is.na(vaf_DF)] = 0
#vaf_DF[vaf_DF >= 0.25] = NA
vaf_mat = apply(X=vaf_DF, MARGIN=2,
                FUN=function(x) as.numeric(x))

# Get columns containing locus info (corresponding to the rows of vaf_mat)
locus_col = full_DF[6:nrow(full_DF), 1]
contig_col = full_DF[6:nrow(full_DF), 2]
start_col = full_DF[6:nrow(full_DF), 3]
end_col = full_DF[6:nrow(full_DF), 4]
gene_col = full_DF[6:nrow(full_DF), 5]
featID_col = full_DF[6:nrow(full_DF), 6]

# Get rows containing sample info (corresponding to the columns of vaf_mat)
sampleID_row = as.character(full_DF[1, 7:ncol(full_DF)])
subjectID_row = as.character(full_DF[2, 7:ncol(full_DF)])
category_row = as.character(full_DF[3, 7:ncol(full_DF)])
run_row = as.character(full_DF[4, 7:ncol(full_DF)])


# Apply partitioning around medoids (PAM) or its extension, clustering for large applications (CLARA),
# rather than k-means clustering because, among other reasons given at the URLs below:
# - from cluster::pam docs: "it is more robust because it minimizes a sum of dissimilarities
# instead of a sum of squared euclidean distances"
# - methods for estimating the optimal number of medoids are more robust than those for estimating
# the optimal number of centroids for data partitioning 
# - medoids are less sensitive to outliers than are centroids obtained by k-means
# - the centroid is an imaginary point in the data (whereas the medoid--the most centrally located data point
# in the cluster--is real), and the initial centroid positions can signficiantly affect the results
# - PAM can use distances other than Euclidean, which is inappropriate for high-dimensional data
# - PAM performs better where expected clusters differ in size and density
# See:
# https://www.kdnuggets.com/2020/12/algorithms-explained-k-means-k-medoids-clustering.html
# https://stats.stackexchange.com/questions/99171/why-is-euclidean-distance-not-a-good-metric-in-high-dimensions
# https://towardsdatascience.com/a-deep-dive-into-partitioning-around-medoids-a77d9b888881

# Estimate the optimal number of clusters for partitioning around medoids (PAM)
# fpc::pamk() criterion="multiasw" with usepam=F (CLARA rather than PAM is used) is
# "recommended for large datasets with 2,000 or more observations; dissimilarity matrices can not be used with clara)"

# Conditionally set pamk_criterion and use_pam
# based on the number of observations in vaf_mat
n_obs = nrow(vaf_mat) * ncol(vaf_mat)
print(n_obs)
if(n_obs < 2000) {
    pamk_criterion = "asw"
    use_pam = T
} else if(n_obs >= 2000) {
    pamk_criterion = "multiasw"
    use_pam = F # CLARA used instead
}

set.seed(4849345)
pamk_n_cl = function(mat, k_range, pamk_criterion, use_pam) {
    fpc::pamk(data=mat,
              krange=k_range,
              criterion=pamk_criterion,
              usepam=use_pam,
              scaling=F,
              alpha=0.001,
              diss=F,
              critout=T,
              ns=2)
}

set.seed(4849345)
vaf_mat_pamk = pamk_n_cl(mat=t(vaf_mat),
                         k_range=1:20,
                         pamk_criterion=pamk_criterion,
                         use_pam=use_pam)
vaf_mat_pamk_n_cl = vaf_mat_pamk$nc
print(paste0("fpc::pamk-estimated optimal number of clusters: ", vaf_mat_pamk_n_cl))
# Use Manhattan rather than Euclidean in view of high dimensionality;
# https://www.kdnuggets.com/2020/12/algorithms-explained-k-means-k-medoids-clustering.html
# https://stats.stackexchange.com/questions/99171/why-is-euclidean-distance-not-a-good-metric-in-high-dimensions
vaf_mat_pam = cluster::pam(x=t(vaf_mat),
                           k=vaf_mat_pamk_n_cl,
                           diss=F,
                           metric="manhattan",
                           stand=F,
                           cluster.only=F,
                           do.swap=T,
                           pamonce=0)

vaf_mat_clara = cluster::clara(x=t(vaf_mat),
                               k=vaf_mat_pamk_n_cl,
                               metric="manhattan",
                               stand=F,
                               cluster.only=F,
                               samples=50,
                               pamLike=T)

# cluster::clara() produces same clusters as cluster::pam()
print(identical(vaf_mat_pam$clustering, vaf_mat_clara$clustering))

pam_row = paste0("Cluster ", vaf_mat_pam$clustering)


# Define colour function for heatmap 
#if(grepl("HIGH", vaf_tsv)) {
#    #quantile_vector = c(0.95, 0.96, 0.97, 0.98, 0.99, 1.00)
#    quantile_vector = c(0.93, 0.94, 0.95, 0.96, 0.97, 0.98)
#} else if(grepl("not_LOW", vaf_tsv)) {
#    quantile_vector = c(0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81)
#}

quantile_vector = c(0.93, 0.94, 0.95, 0.96, 0.97, 0.98)

vaf_mat_colFun_viridis = colorRamp2(quantile(vaf_mat,
                                             quantile_vector,
                                             na.rm=T),
                                    viridis(length(quantile_vector)))
vaf_mat_colFun_plasma = colorRamp2(quantile(vaf_mat,
                                            quantile_vector,
                                            na.rm=T),
                                   plasma(length(quantile_vector)))
vaf_mat_colFun_magma = colorRamp2(quantile(vaf_mat,
                                           quantile_vector,
                                           na.rm=T),
                                  magma(length(quantile_vector)))
# Linear interpolation of colours from VAF=0.00 to VAF=0.25 below obscures VAFs
# at the lower end of the spectrum, hence use of quantile() approach
#vaf_mat_colFun_plasma = colorRamp2(c(0.00, 0.05, 0.10, 0.20, 0.25),
#                                   plasma(5))


# Conditionally set heatmap_size_parameters
# based on the number of rows in vaf_mat
n_rows = nrow(vaf_mat)
if(n_rows < 500) {
    heatmap_width = unit(2, "npc")
    heatmap_height = unit(4, "npc")
    heatmap_height_run = unit(4/nrow(vaf_mat), "npc")
    heatmap_height_category = unit(4/nrow(vaf_mat), "npc")
    pdf_width = 14
    pdf_height = 10
} else if(n_rows >= 500) {
    heatmap_width = unit(40, "cm")
    heatmap_height = unit(160, "cm")
    heatmap_height_run = unit(160/nrow(vaf_mat)*2, "cm")
    heatmap_height_category = unit(160/nrow(vaf_mat)*22, "cm")
    pdf_width = 16
    pdf_height = 68
}


# VAF heatmap
htmp = Heatmap(matrix=vaf_mat,
    col=vaf_mat_colFun_plasma,
    column_split=factor(pam_row, levels=sort(unique(as.character(pam_row)))),
    row_order=1:nrow(vaf_mat),
    row_labels=locus_col,
    row_names_side="left",
    row_names_gp=gpar(fontsize=6),
    show_row_names=TRUE,
    column_title=NULL,
    column_title_rot=0,
    column_title_gp=gpar(fontsize=13),
    column_labels=category_row,
    column_names_rot=90,
    column_names_centered=FALSE,
    cluster_columns=TRUE,
    cluster_column_slices=TRUE,
    clustering_distance_columns="manhattan",
    clustering_method_columns="complete",
    column_dend_side=c("top"),
    #column_dend_height=unit(0.25, "npc"),
    column_dend_height=unit(5, "cm"),
    cluster_rows=TRUE,
    clustering_distance_rows="manhattan",
    clustering_method_rows="complete",
    row_dend_side=c("right"),
    row_dend_width=unit(2, "cm"),
    cluster_row_slices=FALSE,
    heatmap_legend_param=list(title="VAF",
        title_position="topcenter",
        title_gp=gpar(font=2, fontsize=12),
        legend_direction="horizontal",
        labels_gp=gpar(fontsize=10)),
    heatmap_width=heatmap_width,
    heatmap_height=heatmap_height,
    column_gap=unit(1.0, "mm"),
    row_gap=unit(1.0, "mm"),
    row_title=NULL,
    border=FALSE,
    use_raster=TRUE)

## Get column order
#htmp_draw = draw(htmp)
#htmp_column_order = column_order(htmp_draw)


# Define colour code for runs
run_row_sort_uniq = sort(unique(run_row))
run_row_sort_uniq_colours = brewer.pal(n=length(run_row_sort_uniq), "Set2")
names(run_row_sort_uniq_colours) = run_row_sort_uniq

# Run heatmap
run_htmp = Heatmap(mat=matrix(run_row, nrow=1),
    col=run_row_sort_uniq_colours,
    column_split=factor(pam_row, levels=sort(unique(as.character(pam_row)))),
    #column_order=htmp_column_order,
    row_title=NULL,
    row_title_rot=0,
    row_title_gp=gpar(fontsize=13),
    column_labels=NULL,
    cluster_columns=FALSE,
    cluster_column_slices=FALSE,
    cluster_rows=FALSE,
    cluster_row_slices=FALSE,
    heatmap_legend_param=list(title="Run",
        title_position="topcenter",
        title_gp=gpar(font=2, fontsize=12),
        legend_direction="horizontal",
        labels_gp=gpar(fontsize=10)),
    heatmap_width=heatmap_width,
    heatmap_height=heatmap_height_run,
    column_gap=unit(1.0, "mm"),
    row_gap=unit(1.0, "mm"),
    border=FALSE,
    rect_gp=gpar(col="white", lwd=1),
    use_raster=TRUE)


# Category heatmap
category_htmp = Heatmap(mat=matrix(category_row, nrow=1),
    col=c("CRC+"="#be3000", "Polyps+"="#c4c827", "Control"="#008ebe"),
    column_split=factor(pam_row, levels=sort(unique(as.character(pam_row)))),
    #column_order=htmp_column_order,
    row_title=NULL,
    row_title_rot=0,
    row_title_gp=gpar(fontsize=13),
    column_labels=subjectID_row,
    column_names_rot=90,
    column_names_centered=FALSE,
    cluster_columns=FALSE,
    cluster_column_slices=FALSE,
    cluster_rows=FALSE,
    cluster_row_slices=FALSE,
    heatmap_legend_param=list(title="Category",
        title_position="topcenter",
        title_gp=gpar(font=2, fontsize=12),
        legend_direction="horizontal",
        labels_gp=gpar(fontsize=10)),
    heatmap_width=heatmap_width,
    heatmap_height=heatmap_height_category,
    column_gap=unit(1.0, "mm"),
    row_gap=unit(1.0, "mm"),
    border=FALSE,
    rect_gp=gpar(col="white", lwd=1),
    use_raster=TRUE)


htmps = htmp %v% run_htmp %v% category_htmp

legend_gap = unit(15, "mm")

pdf(paste0(plotdir,
           sub(".tsv", "_heatmap.pdf", vaf_tsv)),
    width=pdf_width, height=pdf_height)
draw(htmps,
     gap=unit(1, "mm"),
     column_title=vaf_tsv,
     column_title_gp=gpar(font=1, fontsize=16),
     heatmap_legend_side="bottom",
     legend_gap=legend_gap)
dev.off()
