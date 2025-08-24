#-----------------------------------------------------------------------------------------
# Description: Performs DBSCAN clustering and computes convex hull for cell segmentation 
# Date: Aug 5, 2025
# Author: Hanying Yan
#
# Dependencies:
#   - Seurat
#   - dbscan
#   - dplyr
#   - grDevices 

# Input Data:
#   - Xenium output files
#
# Output Files:
#   - TXT file containing the computed cell boundary in TXT
#-----------------------------------------------------------------------------------------



#######  0. load packages #######  
rm(list=ls())
library(Seurat) 
library(dbscan) 
library(dplyr)
library(grDevices) 


#######  1. load all Xenium samples #######
setwd('path/to/the/project/')
data_dir = "/path/to/the/Xenium/data/"
output_dir = "/path/to/save/outputs/"
gene_name = 'UCHL1' #, 

##########################################################################################
# Function: load_marker
# Description: This function load marker gene from raw Xenium outputs.
# Inputs:
#   - data_dir: the path to where the Xenium data was downloaded
#   - sample_name: name of the specific Xenium sample
#   - gene_name: marker gene name used for cell segmentation, for neuron we use UCHL1
#
# Outputs:
#   - clustering_df: A data frame containing the assigned DBSCAN cluster for each molecule
##########################################################################################
load_marker = function(data_dir, sample_name, gene_name){
  sample_raw = LoadXenium(paste0(data_dir, sample_name))
  sample_input = data.frame(sample_raw@images$fov@molecules$molecules[[gene_name]]@coords)
  return(sample_input)
}

drg1_input = load_marker(data_dir, 'drg1', gene_name) 
drg2_input = load_marker(data_dir, 'drg2', gene_name) 
drg3_input = load_marker(data_dir, 'drg3', gene_name) 
drg4_input = load_marker(data_dir, 'drg4', gene_name) 
tg1_input = load_marker(data_dir, 'tg1', gene_name)



#######  2. define DBSCAN pipeline #######
##########################################################################################
# Function: dbscan_clustering
# Description: This function calculates DBSCAN clustering results given the input data.
# Inputs:
#   - molecules_df: A data frame containing the molecule spatial locations (x, y)
#   - e: selected value for eps
#   - pt: selected value for MinPts
#
# Outputs:
#   - clustering_df: A data frame containing the assigned DBSCAN cluster for each molecule
##########################################################################################
dbscan_clustering = function(molecules_df, e, pt){
  set.seed(1)
  out = dbscan::dbscan(molecules_df, eps = e, minPts = pt)
  clustering_df = molecules_df
  clustering_df$dbscan = as.character(out$cluster)
  n=dim(unique(clustering_df$dbscan))[1]
  print(n)
  return(clustering_df)
}


##########################################################################################
# Function: compute_convex_hull
# Description: This function computes convex hull for the input data.
# Inputs:
#   - clustering_df: A data frame containing the assigned DBSCAN cluster for each molecule
#
# Outputs:
#   - convex_hull_points: A list containing computed convex hull for each DBSCAN cluster
##########################################################################################
compute_convex_hull = function(clustering_df) {
  # compute convex hull given clustering results
  
  # Get cluster IDs and corresponding points
  cluster_ids = clustering_df$dbscan
  cluster_points = clustering_df[,c('x', 'y')]
  
  # Initialize list to store convex hull points
  convex_hull_points = list()
  
  # Iterate over each cluster/cell excluding noise points(dbscan=0)
  for (i in unique(cluster_ids[cluster_ids != 0])) {
    # Get points for current cluster
    cluster_points_i = cluster_points[cluster_ids == i, ]
    # Compute convex hull
    ch = chull(cluster_points_i)
    # Extract convex hull points, adding the first point to complete a polygon
    convex_hull_points_i = cluster_points_i[c(ch, ch[1]),]
    # generate a row starting with 0, the xi,yi 
    convex_hull_row = c(0, as.vector(t(convex_hull_points_i)))
    convex_hull_points[[i]] = convex_hull_row
  }
  return(convex_hull_points)
}


##########################################################################################
# Function: write_convex_hull_to_txt
# Description: This function saves computed convex hulls as TXT.
# Inputs:
#   - sample: sample name
#   - convex_hull_points: A list containing computed convex hull for each DBSCAN cluster
#   - e: selected value for Eps
#   - pt: selected value for MinPts
#   - output_dir: path to save txt file
##########################################################################################
write_convex_hull_to_txt = function(sample, convex_hull_points, e, pt, output_dir) {
  # archive results as txt
  lines = sapply(convex_hull_points, function(row) paste(row, collapse = " "))
  name = paste0(output_dir, sample, '_e', e, '_n', pt, '.txt')
  print(name)
  writeLines(lines, con = name)
}


##########################################################################################
# Function: run_dbscan_save_txt
# Description: This function runs DBSCAN clustering then saves convex hulls as TXT
# Inputs:
#   - molecules_df: A data frame containing the molecule spatial locations (x, y)
#   - e: selected value for Eps
#   - pt: selected value for MinPts
#   - sample: sample name
#   - output_dir: path to save txt file
#
# Outputs:
#   - convex_hull_points: A list containing computed convex hull for each DBSCAN cluster
##########################################################################################
run_dbscan_save_txt = function(molecules_df, e, pt, sample, output_dir) {
  clustering_df = dbscan_clustering(molecules_df, e, pt)
  convex_hull_points = compute_convex_hull(clustering_df)
  write_convex_hull_to_txt(sample, convex_hull_points, e, pt, output_dir)
  print(paste("Finish running:", sample, e, pt, sep=' '))
  return(convex_hull_points)
}



#######  3. test different sets of parameters combination #######
{
  # for all 5 samples: e4n10-30; e5n14-40; e6n24-44, e7n30-56, e8n34-64 (skip = 2)
  # e9n50:120; e10n60:144(skip = 4), e14n80:240(skip = 8), e18n200:520(skip = 16) 
  # total N = 11+14+11+14+19+18+22+21+21 = 151 for drg1-4
  e_vals_all = c(rep(4, 11), rep(5, 14), rep(6, 11), rep(7, 14), rep(8, 19),
                 rep(9, 18), rep(10, 22), rep(14, 21), rep(18, 21))
  pt_vals_all = c(seq(10,30,2), seq(14,40,2), seq(24,44,2), seq(30,56,2), seq(34,70,2),
                  seq(50,120,4), seq(60,144,4), seq(80,240,8), seq(200,520,16))
  
  
  drg1_list = mapply(function(e, pt) {
    run_dbscan_save_txt(drg1_input, e, pt, 'drg1', output_dir)}, 
    e_vals_all, pt_vals_all, SIMPLIFY = FALSE)
  
  drg2_list = mapply(function(e, pt) {
    run_dbscan_save_txt(drg2_input, e, pt, 'drg2', output_dir)}, 
    e_vals_all, pt_vals_all, SIMPLIFY = FALSE)
  
  drg3_list = mapply(function(e, pt) {
    run_dbscan_save_txt(drg3_input, e, pt, 'drg3', output_dir)}, 
    e_vals_all, pt_vals_all, SIMPLIFY = FALSE)
  
  drg4_list = mapply(function(e, pt) {
    run_dbscan_save_txt(drg4_input, e, pt, 'drg4', output_dir)}, 
    e_vals_all, pt_vals_all, SIMPLIFY = FALSE)
  
  tg_list = mapply(function(e, pt) {
    run_dbscan_save_txt(tg1_input, e, pt, 'tg1', output_dir)}, 
    e_vals_all, pt_vals_all, SIMPLIFY = FALSE)
  
  
  # For TG specifically, additional e7n52:60, e8n66:98 (skip = 2)
  # total N = 151+5+14 =170 for tg1
  
  e_vals_tg1 = c(rep(7, 5), rep(8, 14))
  pt_vals_tg1 = c(seq(52, 60, 2), seq(72, 98, 2))
  
  tg_list2 = mapply(function(e, pt) {
    run_dbscan_save_txt(tg1_input, e, pt, 'tg')}, 
    e_vals_tg1, pt_vals_tg1, SIMPLIFY = FALSE)
}


#######  4. decide number of rois to select for new samples and test #######  
##a. load manual annotation for drg1, drg4 and tg1
rois_drg1 = read.csv("trascript_drg_1_roi.csv", row.names = 1) 
rois_drg4 = read.csv("trascript_drg_4_roi.csv", row.names = 1) 
rois_tg1 = read.csv("trascript_tg_1_roi.csv", row.names = 1) 


##########################################################################################
# Function: calculate_range_e
# Description: This function estimates the Eps range to be tested
# Inputs:
#   - rois: A data frame containing estimate size and density of manually annotated cells 
#   - r_median: the median values of estimate cell size
#
# Outputs:
#   - e_range: A vector containing the calculated range of Eps to be tested
##########################################################################################
calculate_range_e = function(rois, r_median){
  lower_e = 0.2*r_median
  upper_e = 0.5*r_median
  e_range = round(lower_e):round(upper_e)
  return (e_range)
}


##########################################################################################
# Function: calculate_range_pt_on_e
# Description: This function estimates the MinPts range to be tested
# Inputs:
#   - e: Eps value to be tested
#   - d_median: the median values of estimate cell density
#
# Outputs:
#   - pt_range: A vector containing the min and max number of MinPts to be tested
##########################################################################################
calculate_range_pt_on_e = function(e, d_median){
  lower_pt = floor(0.6*pi*e*e*d_median)
  upper_pt = ceiling(0.9*pi*e*e*d_median)
  pt_range= c(lower_pt, upper_pt)
  return (pt_range)
}


##########################################################################################
# Function: calculate_range_e_and_pt
# Description: This function estimates the range of Eps and MinPts range to be tested
# Inputs:
#   - rois: A data frame containing estimate size and density of manually annotated cells 
#
# Outputs:
#   - ranges_df: A data frame containing the ranges of Eps and MinPts to be tested
##########################################################################################
calculate_range_e_and_pt = function(rois){
  r_median = median(rois$r)
  d_median = median(rois$density)
  e_range = calculate_range_e(rois, r_median)
  ranges_df = as.data.frame(do.call(rbind, lapply(e_range, function(e) {
    pt_range = calculate_range_pt_on_e(e, d_median)
    c(e, pt_range)
  })))
  
  colnames(ranges_df) = c('Eps', 'lower_MinPts', 'upper_MinPts')
  return(ranges_df)
}


##########################################################################################
# Function: clean_rois
# Description: This function extracts and calculates useful information given rois
# Inputs:
#   - rois: A data frame containing estimate size and density of manually annotated cells 
#
# Outputs:
#   - clean_df: A data frame containing cleaned cell size and density 
##########################################################################################
clean_rois= function(rois){
  #extract rois to get useful information
  clean_df = data.frame(area = rois[,'Area'], n = rois[,'UCHL1'], id = rois[,'roi'], 
                        r = sqrt(rois[,'Area']/pi), r2=rois[,'Area']/pi)
  clean_df$density = clean_df$n/clean_df$area
  return(clean_df)
}


##########################################################################################
# Function: calculate_accuracy
# Description: This function calculates the accuracy of estimated range of Eps and MinPts
# Inputs:
#   - clean_df: A data frame containing extracted information
#   - optimal_e: optimal Eps values 
#   - optimal_pt: optimal MinPts values 
#   - s: number of the selected cells (default = 20)
#   - n_rep: number of replications (default = 10000)
#
# Outputs:
#   - e_and_pt_accuracy: A data frame containing the accuracy of Eps and MinPts per repeat
##########################################################################################
calculate_accuracy = function(rois, optimal_e, optimal_pt, s=20, n_rep=10000){
  
  clean_df = clean_rois(rois)
  
  e_and_pt_accuracy = NULL
  for (seed in 1:n_rep){
    set.seed(seed)
    subsample_df = clean_df[sample(1:nrow(clean_df), s),]
    ranges_df = calculate_range_e_and_pt(subsample_df)
    
    # if optimal Eps is in the range, test whether MinPts MinPts is also in the range
    e_accuracy = ifelse(optimal_e %in% ranges_df$Eps, 'e_in', 'e_out')
    if (e_accuracy == 'e_in'){
      pt_accuracy = ifelse(
        (ranges_df[ranges_df$Eps==optimal_e,'lower_MinPts']<=optimal_pt) & 
          (ranges_df[ranges_df$Eps==optimal_e,'upper_MinPts']>=optimal_pt), 'pt_in', 'pt_out')
    } else (pt_accuracy = "NA")
    e_and_pt = c(e_accuracy, pt_accuracy)
    e_and_pt_accuracy = rbind(e_and_pt_accuracy, e_and_pt)
  }
  
  e_and_pt_accuracy = data.frame(e_and_pt_accuracy)
  colnames(e_and_pt_accuracy) = c("e_range", "pt_range")
  rownames(e_and_pt_accuracy) = paste0(rep("rep", n_rep), 1:n_rep)
  return(e_and_pt_accuracy)
}


##########################################################################################
# Function: test_sample_range
# Description: This function calculates the accuracy given a range of subsampling size
# Inputs:
#   - rois: A data frame containing estimate size and density of manually annotated cells 
#   - optimal_e: optimal Eps values 
#   - optimal_pt: optimal MinPts values 
#   - s_range: a range of selected cells number (default = 20)
#   - n_rep: number of replications (default = 10000)
#
# Outputs:
#   - accuracy_df: A data frame containing the accuracy rates given range of cell numbers
##########################################################################################
calculate_accuracy_range = function(rois, optimal_e, optimal_pt, s_range, n_rep=10000){
  accuracy = c()
  for (s in s_range){
    e_and_pt_acc = calculate_accuracy(rois, optimal_e, optimal_pt, s, n_rep)
    accuracy =  append(accuracy, 
                       table(e_and_pt_acc$e_range, e_and_pt_acc$pt_range)['e_in','pt_in'])
  }
  accuracy_df = data.frame(subsample_size = s_range, accuracy = accuracy)
  accuracy_df$accuracy_rate = accuracy_df$accuracy/100
  return(accuracy_df)
}


## b. test the number of cells selected on drg1
accuracy_df_drg1 = calculate_accuracy_range(rois_drg1, 7, 42, seq(5, 60, 5))


## c. test the formula on drg4
#N=306 for drg4 - randonmly pick N=20; optimal combination is Eps = 7 and MinPts = 48
rois_drg4 = c("roi_17", "roi_186", "roi_183", "roi_105", "roi_131", "roi_113", 
              "roi_48", "roi_108", "roi_99", "roi_13", "roi_215", "roi_178", "roi_78", 
              "roi_276", "roi_208", "roi_75", "roi_9", "roi_187", "roi_249", "roi_141")
clean_df_drg4 = clean_rois(rois_drg4)
summary(clean_df_drg4$r)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.656  13.282  17.756  17.883  21.635  49.375 
clean_df_drg4_select = subset(clean_df_drg4, id %in% rois_drg4)
median(clean_df_drg4_select$r) #21.46
median(clean_df_drg4_select$density) #0.43

accuracy_df_drg4 = calculate_range_e_and_pt(clean_df_drg4_select)
accuracy_df_drg4


###### d. test the formula on tg1 ###### 
#N=2278 for tg1 - randonmly pick N=20; optimal combination is Eps = 6 and MinPts = 40
rois_tg1 = c("roi_801", "roi_534", "roi_1864", "roi_164", "roi_2092", "roi_1725", 
             "roi_1117", "roi_821", "roi_1161", "roi_112", "roi_1013", "roi_1354", 
             "roi_2007", "roi_1484", "roi_906", "roi_381", "roi_2253", "roi_1833", 
             "roi_1275", "roi_1545")
clean_df_tg1 = clean_rois(rois_tg1)
summary(clean_df_tg1$r)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4.146  12.691  15.938  16.207  19.311  39.654 
clean_df_tg1_select = subset(clean_df_tg1, id %in% rois_tg1)
median(clean_df_tg1_select$r) #21.67
median(clean_df_tg1_select$density) #0.46

accuracy_df_tg1 = calculate_range_e_and_pt(clean_df_tg1_select)
accuracy_df_tg1

