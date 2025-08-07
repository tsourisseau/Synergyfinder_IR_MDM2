#synergy analysis of Imed's data on combination of MDM2 inhibitors KTX169 and BI907828 with irradiation
#"RT_Cytation_20250515_Imed_AnalysisTS"

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

setwd("/shared/projects/atacseq/Synergyfinder_IR_MDM2")

#BiocManager::install("synergyfinder")
library("synergyfinder")
library("drc")
library("reshape2")
library("tidyverse")
library("dplyr")
library("tidyr")
library("purrr")
library("furrr")
library("ggplot2")
library("ggforce")
library("grid")
library("vegan")
library("gstat")
library("sp")
library("methods")
library("SpatialExtremes")
library("ggrepel")
library("kriging")
library("plotly")
library("stringr")
library("future")
library("mice")
library("lattice")
library("nleqslv")
library("stats")
library("graphics")
library("grDevices")
library("magrittr")
library("pbapply")
library("metR")



#vignette("synergyfinder")


# Retrieve raw data and shape properly
load("/shared/projects/atacseq/Synergyfinder_IR_MDM2/data/RawData_RT_Cytation_20250515_Imed.RData")

#initialise experiment parameters
dose <- list("0Gy", "2Gy", "4Gy", "6Gy", "8Gy")
ttt <- list("KTX169", "BI907828")
replicates <- 3
Conc_KTX169 <- rep(x = c(0, 1.024E-3, 5.12E-3, 0.0256, 0.128, 0.64, 3.2, 16, 80, 400, 0),
                   times = replicates)
Conc_BI907828 <- rep(x = c(0, 0, 0.00128, 0.0064, 0.032, 0.16, 0.8, 4, 20, 100, 500),
                     times = replicates)
cell_lines <- list("IB111", "T1000", "c93T", "c94T")


#generate average 0nM 0Gy for each cell line (combine 0nM 0Gy for the 2 ttt (KTX and BI))
i_ttt <- NULL
i_dose <- NULL
i_cell_line <- NULL
i_replicate <- NULL
LNum <- 3
for (i_cell_line in c(1:length(cell_lines))) {
  points_0nM_0Gy <- c()
  for (i_replicate in c(1:replicates)) {
    points_0nM_0Gy_N <- unlist(RawData[["0Gy"]][LNum, c(2, 12, 13, 14)])
    points_0nM_0Gy <- c(points_0nM_0Gy, points_0nM_0Gy_N)
    LNum <- LNum + 1
  }
  name_points_0nM_0Gy <- paste0("data_0nM_0Gy_", unlist(cell_lines[i_cell_line]))
  assign(x = name_points_0nM_0Gy, value = points_0nM_0Gy)
}


#generate tables
#initialise iteration variables
i_ttt <- NULL
i_dose <- NULL
i_cell_line <- NULL
i_replicate <- NULL

for (i_ttt in c(1:length(ttt))) {            # for each of the 2 ttt
  #fetch the line from raw data
  if (i_ttt==1) {                      # use either left or right half of the 384-well plate
    Conc_range <- Conc_KTX169          # concentration ranges are different for each inhibitor
    a <- 2              # range of column to fetch data for KTX169
    b <- 12             # range of column to fetch data for KTX169
    c <- 2              # first 0nM point for for KTX169
    d <- 12             # second 0nM point for for KTX169
    matrix_ttt <- as.data.frame(Conc_KTX169, ncol = 1)       # initialise "big" matrix for each ttt
    matrix_ttt <- matrix_ttt[order(matrix_ttt[[1]]), , drop = FALSE]
    colnames(matrix_ttt) <- c("KTX169")
  } else {
    Conc_range <- Conc_BI907828
    a <- 13             # range of column to fetch data for BI907828
    b <- 23             # range of column to fetch data for BI907828
    c <- 13             # first 0nM point for for BI907828
    d <- 14             # second 0nM point for for BI907828
    matrix_ttt <- as.data.frame(Conc_BI907828, ncol = 1)       # initialise "big" matrix for each ttt
    matrix_ttt <- matrix_ttt[order(matrix_ttt[[1]]), , drop = FALSE]
    colnames(matrix_ttt) <- c("BI907828")
  }  # else
  
  rownames(matrix_ttt) <- NULL                        
  
  for (i_dose in c(1:length(dose))){                    # for each dose Gy
    #initialise line number for each new dose (i.e. each new plate)
    LNum <- 3                                  # first line used in 384-well plate
    
    for (i_cell_line in c(1:length(cell_lines))){             # for each cell line
      #re-initialise response vectors
      resp_col_4N_std0Gy <- c()                       # Initialise variable for 3 lines (3 replicates)
      
      
      for (i_replicate in c(1:replicates)){           # 3 replicates
        resp_col_std0Gy <- c(unlist(RawData[[i_dose]][LNum, c(a:b)]))  
        resp_col_std0Gy <- resp_col_std0Gy/mean(get(paste0("data_0nM_0Gy_", unlist(cell_lines[i_cell_line]))))*100              # normalised by ttt=0nM at 0Gy
        
        
        
        #aggregate the 4n into one single vector. this will be the response column
        resp_col_4N_std0Gy <- c(resp_col_4N_std0Gy, resp_col_std0Gy)
        #assign data
        name_cell_lines_ttt_dose_4N_std0Gy <- paste0(unlist(cell_lines[i_cell_line]),      # generate the name of the raw data  for 4 lines
                                                     "_",                                  # eg. IB111_KTX169_0Gy_4N
                                                     unlist(ttt[i_ttt]),
                                                     "_",
                                                     unlist(dose[i_dose]),
                                                     "_4N_std0Gy")
        assign(x = name_cell_lines_ttt_dose_4N_std0Gy, value = resp_col_4N_std0Gy)
        
        #iterate line number
        LNum <- LNum+1
      }  # replicate
      
      #Built data frame for each cell line
      matrix <- as.data.frame(Conc_range, nrow = 1)   # concentration ranges are different for each inhibitor
      matrix <- cbind(matrix, resp_col_4N_std0Gy)
      matrix <- matrix[order(matrix[[1]]), ]
      colnames(matrix) <- c(colnames(matrix)[1], 
                            paste0("col_", 
                                   name_cell_lines_ttt_dose_4N_std0Gy))
      rownames(matrix) <- NULL
      #assign matrix
      name_cell_lines_ttt_dose_4N_matrix_std0Gy <- paste0(name_cell_lines_ttt_dose_4N_std0Gy, 
                                                          "_matrix")
      
      assign(x = name_cell_lines_ttt_dose_4N_matrix_std0Gy, value = matrix)
      
      #built ttt matrix ("big" matrix)
      matrix_ttt <- cbind(matrix_ttt, matrix[[2]])
      colnames(matrix_ttt) <- c(colnames(matrix_ttt)[1:c(length(colnames(matrix_ttt))-1)], 
                                paste0("std0Gy_",
                                       name_cell_lines_ttt_dose_4N_std0Gy))
    }   # cell line
  }     # dose Gy
  # assign "big" matrix_ttt
  name_matrix_ttt_std0Gy <- paste0("matrix_",
                                   unlist(ttt[i_ttt]),
                                   "_std0Gy")
  
  assign(x = name_matrix_ttt_std0Gy, value = matrix_ttt)
  
}       # ttt
# end loop



# KTX169 (conc-column)
######
# KTX169 (conc-column)
conc_c <- c(rep(matrix_KTX169_std0Gy$KTX169, 5))
# Dose Gy (conc-row)
conc_r <- c(rep(c(0, 2, 4, 6, 8), each = 33)) #replicates X 11

# T1000 KTX169
#####
# T1000 KTX169
T1000_KTX_IR <- c(matrix_KTX169_std0Gy$std0Gy_T1000_KTX169_0Gy_4N_std0Gy,
                  matrix_KTX169_std0Gy$std0Gy_T1000_KTX169_2Gy_4N_std0Gy,
                  matrix_KTX169_std0Gy$std0Gy_T1000_KTX169_4Gy_4N_std0Gy,
                  matrix_KTX169_std0Gy$std0Gy_T1000_KTX169_6Gy_4N_std0Gy,
                  matrix_KTX169_std0Gy$std0Gy_T1000_KTX169_8Gy_4N_std0Gy)

#block_id
block_id <- rep(1, length(conc_c))
#DrugRow
DrugRow <- rep("IR", length(conc_c))
#DrugCol
DrugCol <- rep("KTX169", length(conc_c))
#conc-r-unit
conc_r_unit <- rep("Gy", length(conc_c))
#conc-c-unit
conc_c_unit <- rep("nM", length(conc_c))
#Built input data matrix
matrix_T1000_KTX_IR <- data.frame(block_id,
                                    DrugRow,
                                    DrugCol,
                                    conc_r,
                                    conc_c,
                                  T1000_KTX_IR,
                                    conc_r_unit,
                                    conc_c_unit)
colnames(matrix_T1000_KTX_IR) <- c("block_id",
                                     "DrugRow",
                                     "DrugCol",
                                     "conc_r",
                                     "conc_c",
                                     "response",
                                     "conc_r_unit",
                                     "conc_c_unit")



class(matrix_T1000_KTX_IR)
class(matrix_T1000_KTX_IR$response)


#Reshape data
reshape_matrix_T1000_KTX_IR <- ReshapeData(data = matrix_T1000_KTX_IR,
                                             data_type = "viability",
                                             impute = FALSE,
                                             impute_method = NULL,
                                             noise = FALSE, #try without to see if it makes difference
                                             seed = 1)

#calculate synergy scores
syn_reshape_matrix_T1000_KTX_IR <- CalculateSynergy(data = reshape_matrix_T1000_KTX_IR,
                                              method = c("ZIP", "HSA", "Bliss", "Loewe"),
                                              Emin = NA,
                                              Emax = NA,
                                              correct_baseline = "non")


#calculate sensitivity
sensitivity_reshape_matrix_T1000_KTX_IR <- CalculateSensitivity(
  data = reshape_matrix_T1000_KTX_IR,
  correct_baseline = "non"
)


#Make Synergy Scores table
sensitivity_columns <- c("block_id", "drug1", "drug2",
                         "ic50_1", "ic50_2",
                         "ri_1", "ri_2",
                         "css1_ic502", "css2_ic501", "css")
Syn_scores_reshape_matrix_T1000_KTX_IR <- syn_reshape_matrix_T1000_KTX_IR$drug_pairs

write.table(Syn_scores_reshape_matrix_T1000_KTX_IR, file = "Syn_scores_reshape_matrix_T1000_KTX_IR.txt")

#Visualisation
####----------FOR PAPER
#Dose-response curve
#for (i in c(1, 2)){
PlotDoseResponseCurve(
  data = syn_reshape_matrix_T1000_KTX_IR,
  plot_block = 1,
  drug_index = 2,
  plot_new = FALSE,
  record_plot = FALSE,
  plot_subtitle = "T1000 KTX",
  text_size_scale = 1
)
#}

PlotDoseResponseCurve(
  data = syn_reshape_matrix_T1000_KTX_IR,
  plot_block = 1,
  drug_index = 1,
  plot_new = TRUE,
  record_plot = FALSE,
  plot_subtitle = "T1000 IR",
  text_size_scale = 1
)


#Two-drugs combination visualisation
#heatmap
#from survival data: with added statistics = sem 
Plot2DrugHeatmap(
  data = reshape_matrix_T1000_KTX_IR,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "response_origin",
  dynamic = FALSE,
  statistic = "sem",
  summary_statistic = c("mean",  "median"),
  plot_title = "Dose Response Matrix
T1000 KTX169 + IR"
)
# % inhibition: with added statistics = ci (95% confidence interval) or sem
Plot2DrugHeatmap(
  data = reshape_matrix_T1000_KTX_IR,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "response",
  dynamic = FALSE,
  statistic = "sem",
  #  summary_statistic = c("mean",  "median"),
  plot_title = "Inhibition Matrix
T1000 KTX169 + IR")

#Synergy Score Heatmap
Plot2DrugHeatmap(
  data = syn_reshape_matrix_T1000_KTX_IR,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "ZIP_synergy",
  dynamic = FALSE,
  #  summary_statistic = c( "quantile_25", "quantile_75"),
  statistic = "ci",
  plot_title = "ZIP synergy Matrix
T1000 KTX169 + IR"  
)

Plot2DrugHeatmap(
  data = syn_reshape_matrix_T1000_KTX_IR,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "HSA_synergy",
  dynamic = FALSE,
  #  summary_statistic = c( "quantile_25", "quantile_75"),
  statistic = "ci",
  plot_title = "HSA synergy Matrix
T1000 KTX169 + IR"  
)

# 2D contour plot
Plot2DrugContour(data = reshape_matrix_T1000_KTX_IR,
                 plot_block = 1,
                 drugs = c(1, 2),
                 plot_value = "response",
                 dynamic = FALSE,
                 summary_statistic = c("mean", "median"),
                 plot_title = "Inhibition contour plot
T1000 KTX169 + IR")

Plot2DrugContour(data = syn_reshape_matrix_T1000_KTX_IR,
                 plot_block = 1,
                 drugs = c(1, 2),
                 plot_value = "ZIP_synergy",
                 dynamic = FALSE,
                 #                 summary_statistic = c("quantile_25", "quantile_75"))
                 plot_title = "ZIP synergy contour plot
T1000 KTX169 + IR")

# 3D surface plot
Plot2DrugSurface(data = reshape_matrix_T1000_KTX_IR,
                 plot_block = 1,
                 drugs = c(1, 2),
                 plot_value = "response",
                 dynamic = FALSE,
                 summary_statistic = c("mean", "quantile_25", "median", "quantile_75"),
                 plot_title = "Inhibition 3D plot
T1000 KTX169 + IR")

Plot2DrugSurface(data = syn_reshape_matrix_T1000_KTX_IR,
                 plot_block = 1,
                 drugs = c(1, 2),
                 plot_value = "ZIP_synergy",
                 dynamic = FALSE,
                 #                 summary_statistic = c("mean", "quantile_25", "median", "quantile_75"),
                 plot_title = "ZIP synergy 3D plot
T1000 KTX169 + IR",
                 text_size_scale = 2)


#Plotting wrapper: save plts directly into working dir
#e.g. plots all dose response plots
PlotDoseResponse(
  data = reshape_matrix_T1000_KTX_IR,
  block_ids = c(1),
  drugs = c(1,2),
  save_file = TRUE,
  file_type = "png",
  file_name = "smthg3",
  width = 12,
  height = 8
)

#this one returns a ggplot object that can be adjusted with "+ theme()" terms
PlotSynergy(
  data = syn_reshape_matrix_T1000_KTX_IR,
  type = "2D",
  method = "ZIP",
  block_ids = c(1),
  drugs = c(1,2),
  plot_title = "ZIP synergy contour plot
T1000 KTX169 + IR",
  summary_statistic = NULL,
  save_file = TRUE,
  file_type = "png",
  file_name = "IB115_AtrAdct_ZIPsynergy_2Dmap.png" 
)

PlotSynergy(data = syn_reshape_matrix_T1000_KTX_IR,
            type = "3D",
            method = "ZIP",
            block_ids = 1,
            drugs = c(1, 2),
            plot_title = "ZIP synergy 3D plot
T1000 KTX169 + IR",
            summary_statistic = NULL,
            text_size_scale = 1.5,
            save_file = TRUE,
            file_type = "png",
            file_name = "test_wraper_3D_1") + theme(axis.title.y = element_text(angle = 0))




###Plots for paper
#Dose response curves
#size of exported 800x600
PlotDoseResponseCurve(
  data = syn_reshape_matrix_T1000_KTX_IR,
  plot_block = 1,
  drug_index = 2,
  plot_new = FALSE,
  record_plot = FALSE,
  plot_subtitle = "T1000 KTX",
  text_size_scale = 1.5)

#size of exported 800x600
PlotDoseResponseCurve(
  data = syn_reshape_matrix_T1000_KTX_IR,
  plot_block = 1,
  drug_index = 1,
  plot_new = TRUE,
  record_plot = FALSE,
  plot_subtitle = "T1000 IR",
  text_size_scale = 1.5)


# % inhibition: with added statistics = ci (95% confidence interval) or sem
#size of exported 800x600
Plot2DrugHeatmap(
  data = reshape_matrix_T1000_KTX_IR,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "response",
  dynamic = FALSE,
  statistic = "sem",
  #  summary_statistic = c("mean",  "median"),
  plot_title = "Dose Response Matrix
T1000 KTX169 + IR",
  text_label_size_scale = 1)


#Synergy Score Heatmap
#size of exported 800x600
Plot2DrugHeatmap(
  data = syn_reshape_matrix_T1000_KTX_IR,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "ZIP_synergy",
  dynamic = FALSE,
  #  summary_statistic = c( "quantile_25", "quantile_75"),
  statistic = "ci",
  plot_title = "ZIP Synergy Score Heatmap
T1000 KTX169 + IR")

Plot2DrugHeatmap(
  data = syn_reshape_matrix_T1000_KTX_IR,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "HSA_synergy",
  dynamic = FALSE,
  #  summary_statistic = c( "quantile_25", "quantile_75"),
  statistic = "ci",
  plot_title = "HSA Synergy Score Heatmap
T1000 KTX169 + IR")

#Synergy Score 2D-map
#size of exported 800x600
Plot2DrugContour(data = syn_reshape_matrix_T1000_KTX_IR,
                 plot_block = 1,
                 drugs = c(1, 2),
                 plot_value = "ZIP_synergy",
                 dynamic = FALSE,
                 #                 summary_statistic = c("quantile_25", "quantile_75"))
                 plot_title = "ZIP Synergy Score 2D-map
T1000 KTX169 + IR")
#Synergy Score 3D-map
#size of exported 800x600
Plot2DrugSurface(data = syn_reshape_matrix_T1000_KTX_IR,
                 plot_block = 1,
                 drugs = c(1, 2),
                 plot_value = "ZIP_synergy",
                 dynamic = FALSE,
                 #                 summary_statistic = c("mean", "quantile_25", "median", "quantile_75"),
                 plot_title = "ZIP Synergy Score 3D-map
T1000 KTX169 + IR",
                 text_size_scale = 1.5)




# IB111 KTX169
#####
# IB111 KTX169
IB111_KTX_IR <- c(matrix_KTX169_std0Gy$std0Gy_IB111_KTX169_0Gy_4N_std0Gy,
                  matrix_KTX169_std0Gy$std0Gy_IB111_KTX169_2Gy_4N_std0Gy,
                  matrix_KTX169_std0Gy$std0Gy_IB111_KTX169_4Gy_4N_std0Gy,
                  matrix_KTX169_std0Gy$std0Gy_IB111_KTX169_6Gy_4N_std0Gy,
                  matrix_KTX169_std0Gy$std0Gy_IB111_KTX169_8Gy_4N_std0Gy)

#block_id
block_id <- rep(1, length(conc_c))
#DrugRow
DrugRow <- rep("IR", length(conc_c))
#DrugCol
DrugCol <- rep("KTX169", length(conc_c))
#conc-r-unit
conc_r_unit <- rep("Gy", length(conc_c))
#conc-c-unit
conc_c_unit <- rep("nM", length(conc_c))
#Built input data matrix
matrix_IB111_KTX_IR <- data.frame(block_id,
                                  DrugRow,
                                  DrugCol,
                                  conc_r,
                                  conc_c,
                                  IB111_KTX_IR,
                                  conc_r_unit,
                                  conc_c_unit)
colnames(matrix_IB111_KTX_IR) <- c("block_id",
                                   "DrugRow",
                                   "DrugCol",
                                   "conc_r",
                                   "conc_c",
                                   "response",
                                   "conc_r_unit",
                                   "conc_c_unit")



class(matrix_IB111_KTX_IR)
class(matrix_IB111_KTX_IR$response)


#Reshape data
reshape_matrix_IB111_KTX_IR <- ReshapeData(data = matrix_IB111_KTX_IR,
                                           data_type = "viability",
                                           impute = FALSE,
                                           impute_method = NULL,
                                           noise = FALSE, #try without to see if it makes difference
                                           seed = 1)

#calculate synergy scores
syn_reshape_matrix_IB111_KTX_IR <- CalculateSynergy(data = reshape_matrix_IB111_KTX_IR,
                                                    method = c("ZIP", "HSA", "Bliss", "Loewe"),
                                                    Emin = NA,
                                                    Emax = NA,
                                                    correct_baseline = "non")


#calculate sensitivity
sensitivity_reshape_matrix_IB111_KTX_IR <- CalculateSensitivity(
  data = reshape_matrix_IB111_KTX_IR,
  correct_baseline = "non"
)


#Make Synergy Scores table
sensitivity_columns <- c("block_id", "drug1", "drug2",
                         "ic50_1", "ic50_2",
                         "ri_1", "ri_2",
                         "css1_ic502", "css2_ic501", "css")
Syn_scores_reshape_matrix_IB111_KTX_IR <- syn_reshape_matrix_IB111_KTX_IR$drug_pairs

write.table(Syn_scores_reshape_matrix_IB111_KTX_IR, file = "Syn_scores_reshape_matrix_IB111_KTX_IR.txt")

#Visualisation
####----------FOR PAPER
#Dose-response curve
#for (i in c(1, 2)){
PlotDoseResponseCurve(
  data = syn_reshape_matrix_IB111_KTX_IR,
  plot_block = 1,
  drug_index = 2,
  plot_new = FALSE,
  record_plot = FALSE,
  plot_subtitle = "IB111 KTX",
  text_size_scale = 1
)
#}

PlotDoseResponseCurve(
  data = syn_reshape_matrix_IB111_KTX_IR,
  plot_block = 1,
  drug_index = 1,
  plot_new = TRUE,
  record_plot = FALSE,
  plot_subtitle = "IB111 IR",
  text_size_scale = 1
)


#Two-drugs combination visualisation
#heatmap
#from survival data: with added statistics = sem 
Plot2DrugHeatmap(
  data = reshape_matrix_IB111_KTX_IR,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "response_origin",
  dynamic = FALSE,
  statistic = "sem",
  summary_statistic = c("mean",  "median"),
  plot_title = "Dose Response Matrix
IB111 KTX169 + IR"
)
# % inhibition: with added statistics = ci (95% confidence interval) or sem
Plot2DrugHeatmap(
  data = reshape_matrix_IB111_KTX_IR,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "response",
  dynamic = FALSE,
  statistic = "sem",
  #  summary_statistic = c("mean",  "median"),
  plot_title = "Inhibition Matrix
IB111 KTX169 + IR")

#Synergy Score Heatmap
Plot2DrugHeatmap(
  data = syn_reshape_matrix_IB111_KTX_IR,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "ZIP_synergy",
  dynamic = FALSE,
  #  summary_statistic = c( "quantile_25", "quantile_75"),
  statistic = "ci",
  plot_title = "ZIP synergy Matrix
IB111 KTX169 + IR"  
)

Plot2DrugHeatmap(
  data = syn_reshape_matrix_IB111_KTX_IR,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "HSA_synergy",
  dynamic = FALSE,
  #  summary_statistic = c( "quantile_25", "quantile_75"),
  statistic = "ci",
  plot_title = "HSA synergy Matrix
IB111 KTX169 + IR"  
)

# 2D contour plot
Plot2DrugContour(data = reshape_matrix_IB111_KTX_IR,
                 plot_block = 1,
                 drugs = c(1, 2),
                 plot_value = "response",
                 dynamic = FALSE,
                 summary_statistic = c("mean", "median"),
                 plot_title = "Inhibition contour plot
IB111 KTX169 + IR")

Plot2DrugContour(data = syn_reshape_matrix_IB111_KTX_IR,
                 plot_block = 1,
                 drugs = c(1, 2),
                 plot_value = "ZIP_synergy",
                 dynamic = FALSE,
                 #                 summary_statistic = c("quantile_25", "quantile_75"))
                 plot_title = "ZIP synergy contour plot
IB111 KTX169 + IR")

# 3D surface plot
Plot2DrugSurface(data = reshape_matrix_IB111_KTX_IR,
                 plot_block = 1,
                 drugs = c(1, 2),
                 plot_value = "response",
                 dynamic = FALSE,
                 summary_statistic = c("mean", "quantile_25", "median", "quantile_75"),
                 plot_title = "Inhibition 3D plot
IB111 KTX169 + IR")

Plot2DrugSurface(data = syn_reshape_matrix_IB111_KTX_IR,
                 plot_block = 1,
                 drugs = c(1, 2),
                 plot_value = "ZIP_synergy",
                 dynamic = FALSE,
                 #                 summary_statistic = c("mean", "quantile_25", "median", "quantile_75"),
                 plot_title = "ZIP synergy 3D plot
IB111 KTX169 + IR",
                 text_size_scale = 2)


#Plotting wrapper: save plts directly into working dir
#e.g. plots all dose response plots
PlotDoseResponse(
  data = reshape_matrix_IB111_KTX_IR,
  block_ids = c(1),
  drugs = c(1,2),
  save_file = TRUE,
  file_type = "png",
  file_name = "smthg3",
  width = 12,
  height = 8
)

#this one returns a ggplot object that can be adjusted with "+ theme()" terms
PlotSynergy(
  data = syn_reshape_matrix_IB111_KTX_IR,
  type = "2D",
  method = "ZIP",
  block_ids = c(1),
  drugs = c(1,2),
  plot_title = "ZIP synergy contour plot
IB111 KTX169 + IR",
  summary_statistic = NULL,
  save_file = TRUE,
  file_type = "png",
  file_name = "IB115_AtrAdct_ZIPsynergy_2Dmap.png" 
)

PlotSynergy(data = syn_reshape_matrix_IB111_KTX_IR,
            type = "3D",
            method = "ZIP",
            block_ids = 1,
            drugs = c(1, 2),
            plot_title = "ZIP synergy 3D plot
IB111 KTX169 + IR",
            summary_statistic = NULL,
            text_size_scale = 1.5,
            save_file = TRUE,
            file_type = "png",
            file_name = "test_wraper_3D_1") + theme(axis.title.y = element_text(angle = 0))




###Plots for paper
#Dose response curves
#size of exported 800x600
PlotDoseResponseCurve(
  data = syn_reshape_matrix_IB111_KTX_IR,
  plot_block = 1,
  drug_index = 2,
  plot_new = FALSE,
  record_plot = FALSE,
  plot_subtitle = "IB111 KTX",
  text_size_scale = 1.5)

#size of exported 800x600
PlotDoseResponseCurve(
  data = syn_reshape_matrix_IB111_KTX_IR,
  plot_block = 1,
  drug_index = 1,
  plot_new = TRUE,
  record_plot = FALSE,
  plot_subtitle = "IB111 IR",
  text_size_scale = 1.5)


# % inhibition: with added statistics = ci (95% confidence interval) or sem
#size of exported 800x600
Plot2DrugHeatmap(
  data = reshape_matrix_IB111_KTX_IR,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "response",
  dynamic = FALSE,
  statistic = "sem",
  #  summary_statistic = c("mean",  "median"),
  plot_title = "Dose Response Matrix
IB111 KTX169 + IR",
  text_label_size_scale = 1)


#Synergy Score Heatmap
#size of exported 800x600
Plot2DrugHeatmap(
  data = syn_reshape_matrix_IB111_KTX_IR,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "ZIP_synergy",
  dynamic = FALSE,
  #  summary_statistic = c( "quantile_25", "quantile_75"),
  statistic = "ci",
  plot_title = "ZIP Synergy Score Heatmap
IB111 KTX169 + IR")

Plot2DrugHeatmap(
  data = syn_reshape_matrix_IB111_KTX_IR,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "HSA_synergy",
  dynamic = FALSE,
  #  summary_statistic = c( "quantile_25", "quantile_75"),
  statistic = "ci",
  plot_title = "HSA Synergy Score Heatmap
IB111 KTX169 + IR")

#Synergy Score 2D-map
#size of exported 800x600
Plot2DrugContour(data = syn_reshape_matrix_IB111_KTX_IR,
                 plot_block = 1,
                 drugs = c(1, 2),
                 plot_value = "ZIP_synergy",
                 dynamic = FALSE,
                 #                 summary_statistic = c("quantile_25", "quantile_75"))
                 plot_title = "ZIP Synergy Score 2D-map
IB111 KTX169 + IR")
#Synergy Score 3D-map
#size of exported 800x600
Plot2DrugSurface(data = syn_reshape_matrix_IB111_KTX_IR,
                 plot_block = 1,
                 drugs = c(1, 2),
                 plot_value = "ZIP_synergy",
                 dynamic = FALSE,
                 #                 summary_statistic = c("mean", "quantile_25", "median", "quantile_75"),
                 plot_title = "ZIP Synergy Score 3D-map
IB111 KTX169 + IR",
                 text_size_scale = 1.5)










######
# BI907828 (conc-column)
conc_c <- c(rep(matrix_BI907828_std0Gy$BI907828, 5))
# Dose Gy (conc-row)
conc_r <- c(rep(c(0, 2, 4, 6, 8), each = 33)) #replicates X 11


# T1000 BI907828
#####
# T1000 BI907828
T1000_BI_IR <- c(matrix_BI907828_std0Gy$std0Gy_T1000_BI907828_0Gy_4N_std0Gy,
                  matrix_BI907828_std0Gy$std0Gy_T1000_BI907828_2Gy_4N_std0Gy,
                  matrix_BI907828_std0Gy$std0Gy_T1000_BI907828_4Gy_4N_std0Gy,
                  matrix_BI907828_std0Gy$std0Gy_T1000_BI907828_6Gy_4N_std0Gy,
                  matrix_BI907828_std0Gy$std0Gy_T1000_BI907828_8Gy_4N_std0Gy)

#block_id
block_id <- rep(1, length(conc_c))
#DrugRow
DrugRow <- rep("IR", length(conc_c))
#DrugCol
DrugCol <- rep("BI907828", length(conc_c))
#conc-r-unit
conc_r_unit <- rep("Gy", length(conc_c))
#conc-c-unit
conc_c_unit <- rep("nM", length(conc_c))
#Built input data matrix
matrix_T1000_BI_IR <- data.frame(block_id,
                                  DrugRow,
                                  DrugCol,
                                  conc_r,
                                  conc_c,
                                  T1000_BI_IR,
                                  conc_r_unit,
                                  conc_c_unit)
colnames(matrix_T1000_BI_IR) <- c("block_id",
                                   "DrugRow",
                                   "DrugCol",
                                   "conc_r",
                                   "conc_c",
                                   "response",
                                   "conc_r_unit",
                                   "conc_c_unit")



class(matrix_T1000_BI_IR)
class(matrix_T1000_BI_IR$response)


#Reshape data
reshape_matrix_T1000_BI_IR <- ReshapeData(data = matrix_T1000_BI_IR,
                                           data_type = "viability",
                                           impute = FALSE,
                                           impute_method = NULL,
                                           noise = FALSE, #try without to see if it makes difference
                                           seed = 1)

#calculate synergy scores
syn_reshape_matrix_T1000_BI_IR <- CalculateSynergy(data = reshape_matrix_T1000_BI_IR,
                                                    method = c("ZIP", "HSA", "Bliss", "Loewe"),
                                                    Emin = NA,
                                                    Emax = NA,
                                                    correct_baseline = "non")

write.table(Syn_scores_reshape_matrix_T1000_BI_IR, file = "Syn_scores_reshape_matrix_T1000_KTX_BI.txt")

#calculate sensitivity
sensitivity_reshape_matrix_T1000_BI_IR <- CalculateSensitivity(
  data = reshape_matrix_T1000_BI_IR,
  correct_baseline = "non"
)


#Make Synergy Scores table
sensitivity_columns <- c("block_id", "drug1", "drug2",
                         "ic50_1", "ic50_2",
                         "ri_1", "ri_2",
                         "css1_ic502", "css2_ic501", "css")
Syn_scores_reshape_matrix_T1000_BI_IR <- syn_reshape_matrix_T1000_BI_IR$drug_pairs



#Visualisation
####----------FOR PAPER
#Dose-response curve
#for (i in c(1, 2)){
PlotDoseResponseCurve(
  data = syn_reshape_matrix_T1000_BI_IR,
  plot_block = 1,
  drug_index = 2,
  plot_new = FALSE,
  record_plot = FALSE,
  plot_subtitle = "T1000 BI",
  text_size_scale = 1
)
#}

PlotDoseResponseCurve(
  data = syn_reshape_matrix_T1000_BI_IR,
  plot_block = 1,
  drug_index = 1,
  plot_new = TRUE,
  record_plot = FALSE,
  plot_subtitle = "T1000 IR",
  text_size_scale = 1
)


#Two-drugs combination visualisation
#heatmap
#from survival data: with added statistics = sem 
Plot2DrugHeatmap(
  data = reshape_matrix_T1000_BI_IR,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "response_origin",
  dynamic = FALSE,
  statistic = "sem",
  summary_statistic = c("mean",  "median"),
  plot_title = "Dose Response Matrix
T1000 BI907828 + IR"
)
# % inhibition: with added statistics = ci (95% confidence interval) or sem
Plot2DrugHeatmap(
  data = reshape_matrix_T1000_BI_IR,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "response",
  dynamic = FALSE,
  statistic = "sem",
  #  summary_statistic = c("mean",  "median"),
  plot_title = "Inhibition Matrix
T1000 BI907828 + IR")

#Synergy Score Heatmap
Plot2DrugHeatmap(
  data = syn_reshape_matrix_T1000_BI_IR,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "ZIP_synergy",
  dynamic = FALSE,
  #  summary_statistic = c( "quantile_25", "quantile_75"),
  statistic = "ci",
  plot_title = "ZIP synergy Matrix
T1000 BI907828 + IR"  
)

Plot2DrugHeatmap(
  data = syn_reshape_matrix_T1000_BI_IR,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "HSA_synergy",
  dynamic = FALSE,
  #  summary_statistic = c( "quantile_25", "quantile_75"),
  statistic = "ci",
  plot_title = "HSA synergy Matrix
T1000 BI907828 + IR"  
)

# 2D contour plot
Plot2DrugContour(data = reshape_matrix_T1000_BI_IR,
                 plot_block = 1,
                 drugs = c(1, 2),
                 plot_value = "response",
                 dynamic = FALSE,
                 summary_statistic = c("mean", "median"),
                 plot_title = "Inhibition contour plot
T1000 BI907828 + IR")

Plot2DrugContour(data = syn_reshape_matrix_T1000_BI_IR,
                 plot_block = 1,
                 drugs = c(1, 2),
                 plot_value = "ZIP_synergy",
                 dynamic = FALSE,
                 #                 summary_statistic = c("quantile_25", "quantile_75"))
                 plot_title = "ZIP synergy contour plot
T1000 BI907828 + IR")

# 3D surface plot
Plot2DrugSurface(data = reshape_matrix_T1000_BI_IR,
                 plot_block = 1,
                 drugs = c(1, 2),
                 plot_value = "response",
                 dynamic = FALSE,
                 summary_statistic = c("mean", "quantile_25", "median", "quantile_75"),
                 plot_title = "Inhibition 3D plot
T1000 BI907828 + IR")

Plot2DrugSurface(data = syn_reshape_matrix_T1000_BI_IR,
                 plot_block = 1,
                 drugs = c(1, 2),
                 plot_value = "ZIP_synergy",
                 dynamic = FALSE,
                 #                 summary_statistic = c("mean", "quantile_25", "median", "quantile_75"),
                 plot_title = "ZIP synergy 3D plot
T1000 BI907828 + IR",
                 text_size_scale = 2)


#Plotting wrapper: save plts directly into working dir
#e.g. plots all dose response plots
PlotDoseResponse(
  data = reshape_matrix_T1000_BI_IR,
  block_ids = c(1),
  drugs = c(1,2),
  save_file = TRUE,
  file_type = "png",
  file_name = "smthg3",
  width = 12,
  height = 8
)

#this one returns a ggplot object that can be adjusted with "+ theme()" terms
PlotSynergy(
  data = syn_reshape_matrix_T1000_BI_IR,
  type = "2D",
  method = "ZIP",
  block_ids = c(1),
  drugs = c(1,2),
  plot_title = "ZIP synergy contour plot
T1000 BI907828 + IR",
  summary_statistic = NULL,
  save_file = TRUE,
  file_type = "png",
  file_name = "IB115_AtrAdct_ZIPsynergy_2Dmap.png" 
)

PlotSynergy(data = syn_reshape_matrix_T1000_BI_IR,
            type = "3D",
            method = "ZIP",
            block_ids = 1,
            drugs = c(1, 2),
            plot_title = "ZIP synergy 3D plot
T1000 BI907828 + IR",
            summary_statistic = NULL,
            text_size_scale = 1.5,
            save_file = TRUE,
            file_type = "png",
            file_name = "test_wraper_3D_1") + theme(axis.title.y = element_text(angle = 0))




###Plots for paper
#Dose response curves
#size of exported 800x600
PlotDoseResponseCurve(
  data = syn_reshape_matrix_T1000_BI_IR,
  plot_block = 1,
  drug_index = 2,
  plot_new = FALSE,
  record_plot = FALSE,
  plot_subtitle = "T1000 BI",
  text_size_scale = 1.5)

#size of exported 800x600
PlotDoseResponseCurve(
  data = syn_reshape_matrix_T1000_BI_IR,
  plot_block = 1,
  drug_index = 1,
  plot_new = TRUE,
  record_plot = FALSE,
  plot_subtitle = "T1000 IR",
  text_size_scale = 1.5)


# % inhibition: with added statistics = ci (95% confidence interval) or sem
#size of exported 800x600
Plot2DrugHeatmap(
  data = reshape_matrix_T1000_BI_IR,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "response",
  dynamic = FALSE,
  statistic = "sem",
  #  summary_statistic = c("mean",  "median"),
  plot_title = "Dose Response Matrix
T1000 BI907828 + IR",
  text_label_size_scale = 1)


#Synergy Score Heatmap
#size of exported 800x600
Plot2DrugHeatmap(
  data = syn_reshape_matrix_T1000_BI_IR,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "ZIP_synergy",
  dynamic = FALSE,
  #  summary_statistic = c( "quantile_25", "quantile_75"),
  statistic = "ci",
  plot_title = "ZIP Synergy Score Heatmap
T1000 BI907828 + IR")

Plot2DrugHeatmap(
  data = syn_reshape_matrix_T1000_BI_IR,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "HSA_synergy",
  dynamic = FALSE,
  #  summary_statistic = c( "quantile_25", "quantile_75"),
  statistic = "ci",
  plot_title = "HSA Synergy Score Heatmap
T1000 BI907828 + IR")

#Synergy Score 2D-map
#size of exported 800x600
Plot2DrugContour(data = syn_reshape_matrix_T1000_BI_IR,
                 plot_block = 1,
                 drugs = c(1, 2),
                 plot_value = "ZIP_synergy",
                 dynamic = FALSE,
                 #                 summary_statistic = c("quantile_25", "quantile_75"))
                 plot_title = "ZIP Synergy Score 2D-map
T1000 BI907828 + IR")
#Synergy Score 3D-map
#size of exported 800x600
Plot2DrugSurface(data = syn_reshape_matrix_T1000_BI_IR,
                 plot_block = 1,
                 drugs = c(1, 2),
                 plot_value = "ZIP_synergy",
                 dynamic = FALSE,
                 #                 summary_statistic = c("mean", "quantile_25", "median", "quantile_75"),
                 plot_title = "ZIP Synergy Score 3D-map
T1000 BI907828 + IR",
                 text_size_scale = 1.5)


























