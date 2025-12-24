rm(list = ls())
run_AddMulEffects <- function(raw_file_path, harmonised_file_path, outdir, batch_var, formula_str, min_batch_size = 2) {
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
  userdata <- read.csv(raw_file_path, stringsAsFactors = FALSE)
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
  userdata <- read.csv(harmonised_file_path, stringsAsFactors = FALSE)
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
  # addTest() -- test for additive batch effects (use full data)
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
  # Prepare filtered datasets for multiplicative tests
  # (remove tiny batches with count <= min_batch_size)
  ###########################################################
  # Raw
  batch_counts_raw <- as.data.frame(table(raw_data$batch), stringsAsFactors = FALSE)
  colnames(batch_counts_raw) <- c("batch", "N")
  tiny_raw <- batch_counts_raw$batch[batch_counts_raw$N <= min_batch_size]
  if (length(tiny_raw) > 0) {
    message(sprintf("Removing %d tiny batch(es) from multiplicative test on raw data (min_batch_size=%d): %s",
                    length(tiny_raw), min_batch_size, paste(tiny_raw, collapse=", ")))
    raw_for_mult <- raw_data[! raw_data$batch %in% tiny_raw, , drop = FALSE]
  } else {
    raw_for_mult <- raw_data
    message("No tiny batches removed from raw data for multiplicative test.")
  }
  
  # Harmonised
  batch_counts_harm <- as.data.frame(table(harmonised_data$batch), stringsAsFactors = FALSE)
  colnames(batch_counts_harm) <- c("batch", "N")
  tiny_harm <- batch_counts_harm$batch[batch_counts_harm$N <= min_batch_size]
  if (length(tiny_harm) > 0) {
    message(sprintf("Removing %d tiny batch(es) from multiplicative test on harmonised data (min_batch_size=%d): %s",
                    length(tiny_harm), min_batch_size, paste(tiny_harm, collapse=", ")))
    harmonised_for_mult <- harmonised_data[! harmonised_data$batch %in% tiny_harm, , drop = FALSE]
  } else {
    harmonised_for_mult <- harmonised_data
    message("No tiny batches removed from harmonised data for multiplicative test.")
  }
  
  ###########################################################
  # multTest() -- test for multiplicative batch effects
  # (use filtered data where tiny batches removed)
  ###########################################################
  multTestTable <- multTest(idvar='subjectID',
                            batchvar=batch_var,
                            features=featurenames,
                            formula=formula_str,
                            ranef='(1|subjectID)',
                            data=raw_for_mult)
  
  multTestTableCombat <- multTest(idvar='subjectID',
                                  batchvar=batch_var,
                                  features=featurenames_harm,
                                  formula=formula_str,
                                  ranef='(1|subjectID)',
                                  data=harmonised_for_mult)
  
  ###########################################################
  # Save results
  ###########################################################
  write.csv(addTestTable, file = file.path(outdir,"add_raw.csv"), row.names = FALSE)
  write.csv(addTestTableCombat, file = file.path(outdir,"add_harm.csv"), row.names = FALSE)
  write.csv(multTestTable, file = file.path(outdir,"mult_raw.csv"), row.names = FALSE)
  write.csv(multTestTableCombat, file = file.path(outdir,"mult_harm.csv"), row.names = FALSE)
  
  # Also save a small report about removed batches
  report <- list(
    min_batch_size = min_batch_size,
    removed_raw_batches = if(length(tiny_raw)>0) as.character(tiny_raw) else character(0),
    removed_harm_batches = if(length(tiny_harm)>0) as.character(tiny_harm) else character(0),
    batch_counts_raw = batch_counts_raw,
    batch_counts_harm = batch_counts_harm
  )
  saveRDS(report, file = file.path(outdir, "multiplicative_batch_filter_report.rds"))
  
  cat("Processing complete. Results saved to:", outdir, "\n")
  invisible(report)
}
