#-----------------------------------------------------------------------------------------
#Description: Calculates and evaluates proposed parameter range given sub-sample size  
# Date: Aug 5, 2025
# Author: Hanying Yan
#
# Input Data:
#   - transcript_roi.zip: gene expression of manually annotated cells for each sample
#
#-----------------------------------------------------------------------------------------

#######  0. load packages and data #######  
rm(list=ls())
transcripts_drg1 = read.csv("trascript_drg_1_roi.csv", row.names = 1) 
transcripts_drg4 = read.csv("trascript_drg_4_roi.csv", row.names = 1) 
transcripts_tg1 = read.csv("trascript_tg_1_roi.csv", row.names = 1) 


#######  2. evaluates accuracy acorss different subsample sizes  #######  
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
accuracy_df_drg1 = calculate_accuracy_range(transcripts_drg1, 7, 42, seq(5, 60, 5))


## c. test the formula on drg4
#N=306 for drg4 - randonmly pick N=20; optimal combination is Eps = 7 and MinPts = 48
rois_drg4 = c("roi_17", "roi_186", "roi_183", "roi_105", "roi_131", "roi_113", 
              "roi_48", "roi_108", "roi_99", "roi_13", "roi_215", "roi_178", "roi_78", 
              "roi_276", "roi_208", "roi_75", "roi_9", "roi_187", "roi_249", "roi_141")
clean_df_drg4 = clean_rois(transcripts_drg4)
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
clean_df_tg1 = clean_rois(transcripts_tg1)
summary(clean_df_tg1$r)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4.146  12.691  15.938  16.207  19.311  39.654 
clean_df_tg1_select = subset(clean_df_tg1, id %in% rois_tg1)
median(clean_df_tg1_select$r) #21.67
median(clean_df_tg1_select$density) #0.46

accuracy_df_tg1 = calculate_range_e_and_pt(clean_df_tg1_select)
accuracy_df_tg1

