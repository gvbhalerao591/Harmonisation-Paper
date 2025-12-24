setwd("/Users/psyc1586_admin/GVB_data/harmonisation_work/00_aa_skullImprovement_Oct2025/run_harmonisation_pipelines/longComBat")
filename = "data_for_longComBat.csv"

# Load your data
data <- read.csv(filename, stringsAsFactors = FALSE)
data <- data[data$subjectID != "UCL006", ] 

# Define IDPs
idp_names <- c("T1_SIENAX_periphGM_norm_vol", 
               "T1_SIENAX_CSF_norm_vol",
               "T1_SIENAX_GM_norm_vol",
               "T1_SIENAX_WM_norm_vol",
               "T1_SIENAX_brain_norm_vol",
               "T1_FIRST_left_hippocampus",
               "T1_FIRST_right_hippocampus")

# Run harmonisation: 1) scanner as fixed effect and site as random effect
result_scannersite <- harmonise_longComBat(
  data = data,
  batch_variable = "scanner_site",  # Scanner fixed + site random
  idp_names = idp_names,
  save_output = TRUE
)

# Run harmonisation: 1) scanner as fixed effect 
result_scanner <- harmonise_longComBat(
  data = data,
  batch_variable = "scanner",  # Scanner fixed + site random
  idp_names = idp_names,
  save_output = TRUE
)

# Run harmonisation: 1) site as fixed effect 
result_site <- harmonise_longComBat(
  data = data,
  batch_variable = "site",  # Scanner fixed + site random
  idp_names = idp_names,
  save_output = TRUE
)
