setwd("/Users/psyc1586_admin/GVB_data/harmonisation_work/00_aa_skullImprovement_Oct2025/run_harmonisation_pipelines/CovBat_approach")

batch_col = "Site"
# Read your data
df <- read.csv("data_forCovBat.csv")
df <- df[df$subjectID != "UCL006", ]  # Remove outlier if needed

# Define your features
feature_cols <- c("T1_SIENAX_periphGM_norm_vol", "T1_SIENAX_CSF_norm_vol", 
                  "T1_SIENAX_GM_norm_vol", "T1_SIENAX_WM_norm_vol", 
                  "T1_SIENAX_brain_norm_vol", "T1_FIRST_left_hippocampus",
                  "T1_FIRST_right_hippocampus")

# Run harmonization
result <- harmonize_CovBat(
  data = df,
  feature_cols = feature_cols,
  batch_col = batch_col,
  age_col = "Final_Age",
  timepoint_col = "timepoint",
  reference_timepoint = "TP0"
)

# Access harmonized data
harmonized_df <- result$harmonized_df

# Save results
write.csv(result$harmonized_df, 
          paste0("covbat_harmonised_", batch_col, "_output.csv"), 
          row.names = FALSE)
# Access full CovBat result if needed for diagnostics
covbat_stats <- result$covbat_result