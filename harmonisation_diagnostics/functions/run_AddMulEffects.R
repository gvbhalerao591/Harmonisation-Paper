rm(list = ls())
run_AddMulEffects <- function(raw_file_path, harmonised_file_path, outdir, batch_var, formula_str) {
  library(longCombat)
  library(invgamma)
  library(lme4)
  
  # Create output directory if it doesn't exist
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  ###########################################################
  # Load Raw Data
  ###########################################################
  userdata <- read.csv(raw_file_path)
  required_cols <- c("subjectID", "timepoint", "zscore_age", batch_var)
  missing_cols <- setdiff(required_cols, names(userdata))
  if (length(missing_cols) > 0) {
    warning(paste("Missing columns in raw data:", paste(missing_cols, collapse=", ")))
  }
  raw_data <- userdata
  featurenames <- grep("^T1", colnames(raw_data), value=TRUE)
  
  ###########################################################
  # Load Harmonised Data
  ###########################################################
  userdata <- read.csv(harmonised_file_path)
  required_cols <- c("subjectID", "timepoint", "zscore_age", batch_var)
  missing_cols <- setdiff(required_cols, names(userdata))
  if (length(missing_cols) > 0) {
    warning(paste("Missing columns in harmonised data:", paste(missing_cols, collapse=", ")))
  }
  harmonised_data <- userdata
  featurenames_harm <- grep("^harmonised", colnames(harmonised_data), value=TRUE)
  
  ###########################################################
  # Print user-provided formula
  ###########################################################
  print(paste("Using formula:", formula_str))
  
  raw_data$batch <- factor(raw_data[[batch_var]])
  harmonised_data$batch <- factor(harmonised_data[[batch_var]])
  ###########################################################
  # addTest() -- test for additive batch effects
  ###########################################################
  addTestTable <- addTest(idvar='subjectID',
                          batchvar=batch_var,
                          features=featurenames,
                          formula=formula_str,
                          ranef='(1|subjectID)',
                          data=raw_data)
  
  addTestTableCombat <- addTest(idvar='subjectID',
                                batchvar=batch_var,
                                features=featurenames_harm,
                                formula=formula_str,
                                ranef='(1|subjectID)',
                                data=harmonised_data)
  
  ###########################################################
  # multTest() -- test for multiplicative batch effects
  ###########################################################
  multTestTable <- multTest(idvar='subjectID',
                            batchvar=batch_var,
                            features=featurenames,
                            formula=formula_str,
                            ranef='(1|subjectID)',
                            data=raw_data)
  
  multTestTableCombat <- multTest(idvar='subjectID',
                                  batchvar=batch_var,
                                  features=featurenames_harm,
                                  formula=formula_str,
                                  ranef='(1|subjectID)',
                                  data=harmonised_data)
  
  ###########################################################
  # Save results
  ###########################################################
  write.csv(addTestTable, file = file.path(outdir,"add_raw.csv"), row.names = FALSE)
  write.csv(addTestTableCombat, file = file.path(outdir,"add_harm.csv"), row.names = FALSE)
  write.csv(multTestTable, file = file.path(outdir,"mult_raw.csv"), row.names = FALSE)
  write.csv(multTestTableCombat, file = file.path(outdir,"mult_harm.csv"), row.names = FALSE)
  
  cat("Processing complete. Results saved to:", outdir, "\n")
}
