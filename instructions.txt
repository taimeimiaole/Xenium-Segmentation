README


Creating the data in order to train the different segmentation softwares is a very demanding task so we have uploaded both the raw data used to train the different software in addition to the prepared data after cleanup. Furthermore, the code used to clean the data is also included. Regarding the models, the files for each Baysor run is provided in addition to code involved in training and testing the Cellpose, YoloV8, DBScan, and Combo models. Code included for further analysis from Xenium is also included in the Xenium_manuscript_code_HY folder. 


Baysor and Training


aligned_he:
Folder which stores all of the aligned H&E images that line up precisely with the ST images generated


alignment_files:
Set of alignment files for the different tissues to transform the pre-aligned H&E images to be aligned with the ST images. Pipeline uses the alignment matrix to transform pre-aligned H&E image to aligned H&E image for the purpose of this project.


Batch2:
Set of csv files which store the manual masks/annotations for each corresponding 500 by 500 pixel pieces of the tissue (for the DRG and TG tissues). Also contains the corresponding pieces of the aligned H&E image and the original ST image that match the csv files. Dataset is used for training both the YoloV8, Cellpose, and combination models used in this project.


Baysor_yaml_files:
This folder contains the data involved in each Baysor run for all the DRG and TG tissues. It involves both the min_molec_per_gene constraints and the scale constraints. To obtain the yaml file used to initiate each Baysor run, each set is stored as the “config.toml” file in the respective trial folders under the larger “min_molec” and “scale” folders for each tissue. 
* The input to the terminal is given in the run_parameters_and_record.odt file within the folder as well. The file paths need to be adjusted based on what dataset is used but all that is needed is an input folder containing a config.toml file as well as an appropriate output folder for all the data (which can be the same as the input folder). Also needed is a set of filtered transcripts for the respective tissue being segmented. 






Crop_Images:
This folder contains the pipeline to crop images for any specific purpose; input is any image along with the top-left and bottom-right pixel coordinates of the image following the cropping. Output is the newly cropped image


csv_to_roi:
This folder contains a pipeline which, given a csv file containing all the manual annotations/masks in the specified format, converts them to a list of rois which are stored in a separate folder. Also contains the code to zip the folder of rois for upload into ImageJ for visualization. Included in the folder is a subfolder with all the csv data along with the unzipped and zipped rois for the manual annotations belonging to each tissue.


Filtered_Xenium_Transcripts:
This folder stores the processed Xenium transcripts obtained by filtering the raw Xenium transcripts for each tissue. After processing, the transcripts can be fed into Baysor for segmentation of the different tissues and corresponding metrics.


Format_Conversion:
This folder contains a set of code that can convert between different formats of file/document that were used throughout this project. Given one specific format and the file/folder path, this pipeline is able to transform it into the desired output/format.


He_stain:
Contains the set of raw H&E images pre-alignment as well as a pipeline to convert the raw H&E image into a properly aligned H&E image. Input required is the original H&E image as well as the corresponding transformation matrix, which was downloaded from Xenium in our project. 


manual_annotations:
Contains the manual_annotations in csv format which were hand-drawn in ImageJ for the different tissues used in this project.


RoiSet:
This folder stores the raw .roi files for each of the tissues used in this project within separate subfolders. Also contains a conversion code which can convert a folder of rois for each tissue into a csv storing the annotations for each tissue that can be fed into other pipelines. The input given should only include the folder storing all the rois, and the output is the corresponding csv file stored at the specified path.








Split_Batch2:
Set of Batch2 training data split into Training, Testing and Validation subfolders (80-10-10) for corresponding H&E and ST images as well as the csv files. 


Split_Data:
Set of different pipelines to split up data for the corresponding 5 tissues used in this study. Input to this pipeline is the Batch2 folder and the output is a new folder containing the original data split randomly into training, testing, and validation folders.


Training_Data_Generation:
Set of different pipelines to generate Batch2 folder; note that each of the 5 pipelines generates the corresponding training set for a separate tissue, so all 5 need to be run to obtain the Batch2 dataset. Input is the corresponding csv file for the manual annotations belonging to each tissue as well as necessary aligned H&E and ST images for the respective tissue. Output is the corresponding subfolders within Batch2 for each of the 5 tissues. 


transcripts:
This folder contains the set which includes all the different ST images for each of the 4 DRG tissues and the singular TG tissue used in this project. Within the ST images, the UCHL1 and FABP7 transcripts are highlighted to show the outline of the tissue.


Xenium_Data:
This folder contains a pipeline to create the ST images stored in the transcripts folder. The input is a folder which contains all of the required Xenium Data and the output is an ST image for the corresponding tissues for the input folder. The input folder containing all the Xenium raw data can be downloaded from the public dataset once the manuscript is accepted for publication. 


YoloV8 and Code


Merged Image Generation: 
This code creates the merged image dataset by first splitting ST and HE images and then merging their respective single channel images. 


YOLO Training: 
Code used for YOLO training. Parameters set include task (detection, segmentation, pose, or classification), imgsz (image size), model size, and epoch number. Input is the DRG123_batch2 dataset for ST model, DRG123_HE_batch2 for HE model, DRG123_Merged for Merged model.




YOLO Prediction: 
Code for using YOLO for prediction. Image size must be set, as well as max_det if the amount of detections in the image is estimated to be over 300. Save_txt is also set to true to save the predictions as output .txt files. Input images may have to be cropped and then stitched back together if the max_det is too large or image size too large.


Cropping Images for Prediction: 
Crops images into specified coordinates for use in the YOLO prediction code.


Stitch: 
Stitches output .txt files of YOLO predictions back together for a singular image, using the same specified coordinates used to crop the image. Input includes the cropped images, original full image, output .txt files for cropped images, and specified output path for stitched .txt file and image.


Training Different Combinations: 
Code used for the initial training of different combinations of architectures and encoders to compare performance and accuracy. Combinations used are specified in combination.txt file with architecture listed first, then encoder, then encoder weights. Input should be separated into train, test, and val (not valid!) folders with image and label subfolders for each. The dataset does not need a dataset.yaml file. Output includes csv logs and tb logs for each version trained. Csv logs include a csv file of val_map50 scores after each training epoch. 


Extracting val_map50 values for Combination Model: 
Extracts all the val_map50 scores for each combination and organizes them in a csv file. Using this csv file in Excel, average val_map50s can be determined and models can be quantitatively compared. All paths to each csv file for each combination can be specified in the code so it only has to be run once.


Combination Model - Binary Segmentation: 
Used for the top 3 combinations for model training and prediction. Prediction is done on images within the test subfolder and must be the same dimensions as training images (specified image size in code). The input dataset is similar to the YOLO input dataset, except the validation folder is specified as “val.” Also, no dataset.yaml file is needed. Output is binary_pr npy arrays and pr npy arrays.   


Combining Combo Model Output Arrays: 
Organizes Combo model output from binary segmentation code into binary_pr npy arrays and pr npy arrays. Binary_pr npy arrays are resized and then combined to form one array for the original full input image.
CellPose


The Cellpose training pipeline is a very extensive process with many parameter specificities that the user must keep track of. We have created a broad set of training, visualization, and evaluation tools to allow for more streamlined use of Cellpose. Find out how to use our code below.


CLEANING


We have given an automated way to compile the datasets necessary for Cellpose training given the cropped DRG and TG folders we have available. The process of cleaning data for usage in Cellpose involves: 


merge_folders_final.ipnyb
* Merging our original dataset folders to compile DRG1-3 data into a unified dataset
* Inputs: List of directories to merge, output directory, handling duplicate data
* Outputs: Merged directory


yolo_csv_conversion_final.ipnyb
* Converting yolo_txt ground-truth files into csv files for use in downstream model analysis
* Inputs: Image size, directory
* Outputs: dir + ‘/csv/’ subfolder of all of the generated csv files


merge_folders_final.ipnyb
* Converting csv ground-truth files into numpy arrays (.pngs) for use in Cellpose training
* Inputs: csv directory generated from the previous file, output (image) directory
* Outputs: ‘_mask.png’ files ready for training


TRAINING


Training in Cellpose has been divided into two steps: Pretrained model selection and parameter optimization (channel selection, n_epochs, min_train_masks, and flow_prob_threshold).


Important Notes:


Our evaluation pipeline is incorporated in every step of training so the path to the evaluation_pipeline code should be adjusted in the %run line at the top of each pipeline.




You may notice an img_folder parameter at the top of each file. This parameter is used to more efficiently switch between HE images and transcriptomics image usages and is based off of the established folder names.
There is an if_else statement that you can modify to use the img_folder parameter. Otherwise, remove the img_folder parameter and modify the parameters in the if_else statement to suit your needs:
* eval_dir: The directory where the output csvs will be for evaluation purposes
* pred_path: The output file that will be evaluated (will be the same name as the input image filename but .csv)
* out_dir: The directory you want the evaluation output to go in
* man_path: The corresponding csv file generated for the particular evaluation image
* For training there is additional information to add:
   * out_path_drg4: Provides the out_path for the evaluation of the model performance on the established evaluation image (DRG4 for our purposes)
   * out_path_test: Provides the out_path for the evaluation of the model performance on the cropped test image folder
We used DRG4 as our evaluation image but this can be changed based on your needs.


Pretrained Model Selection:


Cellpose has many pretrained models available for selection in the models.MODEL_NAMES list. The transformer and neurips models were not tested because it requires a particular back-end set up that is incompatible with all of the other models and were irrelevant to the goal of this step in training. 
If the user has their own pretrained model they want to optimize via training, model_list can be adjusted to incorporate this model.


run_cellpose_default_models_final.ipynb


* Tests all provided pretrained models using default parameters on a given input image
* Inputs:
   * Change img_folder and corresponding if-else statement as described in Important Notes
   * model_list: defaults to models.MODEL_NAMES without neurips or transformer models
   * c1, c2: Red, Blue, Green, Grayscale (See Channel Selection for more information)
* Outputs:
   * Csv folder in the output directory with all of the evaluation outputs
   * Heatmaps to visualize probability maps in the output directory
   * (See EVALUATION for default model comparison and visualization)




run_cellpose_default_training_final.ipynb


* Tests selected pretrained models for further training and selection
   * While it is ideal that the initial default_models analysis can yield the best model for training, it is not always the case that it is the best model. Instead, the goal is to eliminate very poor performing models to see how adaptable these better performing models are to training in this step.
* Inputs:
   * Change img_folder and corresponding if-else statement as described in Important Notes (inside the for loop)
   * initial_models: hard coded as the specific models that performed the best for our DRG data
      * Can be adjusted for user’s purposes
   * base_dir: DRG123 folder (Crucial that inner folders are formatted as are to allow for downstream train and test folders to be identified correctly)
   * txt_folder: The folder that contains the ground truth masks (Cleaned in earlier step)
   * c1, c2: Red, Blue, Green, Grayscale (See Channel Selection for more information)
* Outputs:
   * Trained model files in train_dir under a generated models folder
   * Csv folder in the evaluation directory with evaluation results
   * Heatmaps in heatmap directory to visualize probability maps
   * (See EVALUATION for default model comparison and visualization)


Parameter Evaluation


run_cellpose_channels_final.ipynb
* Identifies the best channel to use for prediction
* It is recommended that this step be performed before and after default channel selection (use cyto2 or cyto3 as the default model) to ensure that channel selection does not impact model selection and vice versa
* Inputs:
   * Change img_folder and corresponding if-else statement as described in Important Notes (inside the for loop)
   * initial_model: Best performing model from pretrained model selection
   * base_dir: DRG123 folder (Crucial that inner folders are formatted as are to allow for downstream train and test folders to be identified correctly)
   * txt_folder: The folder that contains the ground truth masks (Cleaned in earlier step)
   * c1: Red given our neuronal marker & H&E stain were both red
   * c2_list: List of second channel choices (Green, None for transcript images given our glial cell marker was green)
      * Can be adopted to the user’s purpose (c1 can also be changed to be the variable)
* Outputs:
   * Trained model files in train_dir under a generated models folder
   * Csv folder in the evaluation directory with evaluation results
   * Heatmaps in heatmap directory to visualize probability maps
   * (See EVALUATION for channel comparison and visualization)


run_cellpose_epochs_final.ipynb
* Identifies the best epoch choice for prediction
* Inputs:
   * Change img_folder and corresponding if-else statement as described in Important Notes (inside the for loop)
   * initial_model: Best performing model from pretrained model selection
   * base_dir: DRG123 folder (Crucial that inner folders are formatted as are to allow for downstream train and test folders to be identified correctly)
   * txt_folder: The folder that contains the ground truth masks (Cleaned in earlier step)
   * heatmap_dir: Directory of heatmaps (not as significant for this step)
   * c1, c2: Optimized channel selection for the given images
   * e_list: List of epochs to test currently from 1-300
      * Can be adopted to user’s purposes (This is a long step and certain epoch choices can be removed if need be)
* Outputs:
   * Trained model files in train_dir under a generated models folder
   * Csv folder in the evaluation directory with evaluation results
   * Heatmaps in heatmap directory to visualize probability maps
   * (See EVALUATION for epoch comparison and visualization)


run_cellpose_min_masks_final.ipynb
* Identifies the best min_train_masks choice for prediction
* Inputs:
   * Change img_folder and corresponding if-else statement as described in Important Notes (inside the for loop)
   * initial_model: Best performing model from pretrained model selection
   * base_dir: DRG123 folder (Crucial that inner folders are formatted as are to allow for downstream train and test folders to be identified correctly)
   * txt_folder: The folder that contains the ground truth masks (Cleaned in earlier step)
   * heatmap_dir: Directory of heatmaps (not as significant for this step)
   * c1, c2: Optimized channel selection for the given images
   * epochs: Optimized epoch selection for the given images
   * max_masks: Upper bound for minimum mask threshold
      * We currently have it as larger than the maximum number of masks in any given cropped image but this can be adjusted based on user needs
* Outputs:
   * Trained model files in train_dir under a generated models folder
   * Csv folder in the evaluation directory with evaluation results
   * Heatmaps in heatmap directory to visualize probability maps
   * (See EVALUATION for min_mask comparison and visualization)


run_cellpose_flow-prob_final.ipynb
* Identifies the best flow threshold and probability threshold combination for prediction
* Inputs:
   * Change img_folder and corresponding if-else statement as described in Important Notes (inside the for loop)
   * base_dir: DRG123 folder (Crucial that inner folders are formatted as are to allow for downstream train and test folders to be identified correctly)
   * txt_folder: The folder that contains the ground truth masks (Cleaned in earlier step)
   * c1, c2: Optimized channel selection for the given images
   * min_masks: Optimized min_train_mask threshold for the given images
   * model_name: Automatically identified given the previous parameters provided (should exist given the model_name formatting is preserved from the min_mask step)
   * flow_threshold and cellprob_threshold have predefined ranges and step sizes that will be tested
* Outputs:
   * Three csvs of evaluation parameters given every tested combination of flow_prob_thresholds
   * (See EVALUATION for flow-prob comparison and visualization)


run_cellpose_pretrained_final.ipynb
* Predicts masks of a given directory of images
* Inputs:
   * model_name: Pretrained/optimized model to perform the prediction
   * dir: Directory with desired images
   * flow_threshold and cellprob_threshold: Default or optimized values (change in the model.eval call directly)
* Outputs:
   * Predicted ROIs of the images
   * OPTIONAL: Output heatmap


evaluationPipeline.ipynb
* Evaluates the prediction performance based on our established parameters
* Folder based evaluation (optimized_model_eval)
   * Inputs:
      * man_dir: Manual directory (csvs)
      * pred_dir: Prediction directory (csvs)
      * out_path: Path to output csv
   * Outputs:
      * Evaluation csv (See paper for explanation of parameters)
* File based evaluation (optimized_model_eval_files)
   * Inputs:
      * man_path: Ground truth mask path (csv)
      * pred_path: Predicted mask path (csv)
      * Out_path: Path to output csv
   * Outputs:
      * Evaluation csv (See paper for explanation of parameters)
* Notes:
   * This pipeline is tied to the Cellpose training pathway and is run directly in training so its functionality is currently commented out to allow for streamlined Cellpose training
   * Follow instructions in the pipeline itself to use it directly


EVALUATION


After every step in training, the results should be evaluated to allow for more optimized training in the subsequent steps. There are three main ways that this evaluation occurs: Probability maps (direct output of training), stacked bar graphs (can be relative), and F1-score scatter plots


Probability Maps: 
Qualitative way to evaluate the model performance to see if the probability generally correlates with the expected cell locations


Stacked Bar Graphs: 
Mixed method to quantify the specific type of ground-truth-prediction interactions (i.e. true positives, false negatives, partial overlap, encompassing prediction, etc.)
These bar graphs are based on the total number of interactions between ground truth masks and predicted masks so the total number of interactions for each model may not be the same. Consequently, we also normalize the bar graphs to better compare model performances based on percentage of each categorization rather than just the number of such categorizations.


We have two types of bar graph sorting methods:
value_sorted_bar_final.ipynb
* This is used for default model comparison, default model training comparison, and channel comparison (If >2 channel combinations are used)
* Sorts based on the total number of true positives predicted
* Inputs:
   * name: Name of the figure
   * comparison_folder: Input (and output) folder with all of the evaluation csvs
* Outputs:
   * Stacked bar graph and normalized bar graph
name_sorted_bar_final.ipynb
* This is used for DRG image comparison (evaluation of model performance on every full DRG image), model comparison (between Cellpose, Baysor, and YOLOv8), and channel comparison (If <=2 channel combinations are used)
* Sorts based on the name of the file
* Inputs:
   * name: Name of the figure
   * comparison_folder: Input (and output) folder with all of the evaluation csvs
* Outputs:
   * Stacked bar graph and normalized bar graph


F1-Score Scatter Plots: 
Quantitative method to track the f1-scores across a range of values and identify potential trends and peaks.


scatter_plot_final.ipynb
* This is used for epoch and min_mask evaluation and both may have different parameters but overall have the same set up:
* Inputs:
   * name: Name of the figure
   * batch: Batch 2 was used for all evaluation
   * ev_imgs: Can be drg4 or test based on which the user is more interested in
   * dir: Directory to the evaluation csvs


Flow-Prob Evaluation: We made a 3D-plot to visualize the metrics of each flow-prob combination and were particularly interested in F1-score.
flow-prob_3D-plot_final.ipynb
* Evaluates the prediction performance based on our established parameters
* Folder based evaluation (optimized_model_eval)
   * Inputs:
      * model_name: The type of model used (for output path naming purposes)
      * input_path: Path to F1-50 csv (can also use the F1-50-95 csv)
   * Outputs:
      * 3D-plot of F1-50 scores with different orientations for different visualization


Xenium_manuscript_code_HY


This folder contains three subfolders, which contain: one R Markdown file and two Jupyter Notebook files. These scripts are used for downstream analyses, including cell clustering, cell type correlation analysis, marker gene heatmap generation, and dot plot visualization.


generate_gene_expression_matrix_from_model_segmentation: This notebook generates a gene expression matrix based on the model segmentation output.


co-clustering_and_cell_type_correlation_analysis: This R Markdown file uses Seurat to perform cell clustering and related analyses.


xenium_marker_gene_plot: This notebook generates dot plots to visualize marker gene expression.


Xenium_manuscript_code_DBSCAN


This folder contains two subfolders:


1. DBSCAN_cell_segmentation:
DBSCAN_cell_segmentation.R: This R script loads raw Xenium outputs (submitted to the Pennsieve dataset) to extract marker gene expression. It then applies DBSCAN clustering and computes the convex hull for each cluster to define cell boundaries. The resulting cell segmentation is saved as a .txt file, which can be used in our evaluation pipeline for accuracy assessment.


2. Proposed_parameter_range_accuracy:
a). Proposed_parameter_range_accuracy.R: This R script evaluates whether the parameter range defined by our proposed formula includes the optimal combination of Eps and MinPts, based on iterative subsampling tests, across different subsample sizes. It also applies the formula to new samples and assesses whether it captures appropriate parameter settings.
b). transcript_roi.zip: The input data for Proposed_parameter_range_accuracy.R, including drg_1, drg_4 and tg_1.
