harmonize_CovBat <- function(data, 
                                  feature_cols, 
                                  batch_col,
                                  age_col,
                                  timepoint_col,
                                  reference_timepoint = "TP0",
                                  metadata_cols = c("subjectID", "timepoint", "age", "Site", "Scanner", "groupname"),
                                  percent_var = 0.95,
                                  standardize_var = TRUE,
                                  mean_only = FALSE,
                                  empirical_bayes = TRUE,
                                  parametric = TRUE,
                                  verbose = TRUE) {
  
  # Load required library
  if (!require(CovBat)) {
    stop("CovBat package is required. Please install it first.")
  }
  
  # Validate inputs
  if (!all(feature_cols %in% colnames(data))) {
    stop("Not all feature columns are present in the data")
  }
  if (!(batch_col %in% colnames(data))) {
    stop("Batch column not found in data")
  }
  if (!(age_col %in% colnames(data))) {
    stop("Age column not found in data")
  }
  if (!(timepoint_col %in% colnames(data))) {
    stop("Timepoint column not found in data")
  }
  
  # Create a copy to avoid modifying original data
  df <- data
  
  # Create data matrix (features x subjects)
  dat <- t(as.matrix(df[, feature_cols]))
  
  # Z-score age
  df$age <- scale(df[[age_col]])[, 1]
  
  # Create batch variable
  bat <- as.factor(df[[batch_col]])
  
  # Create model matrix of covariates (no intercept)
  mod_full <- model.matrix(~ 0 + age + timepoint, data = df)
  
  # Drop reference timepoint dummy column
  reference_col <- paste0(timepoint_col, reference_timepoint)
  if (reference_col %in% colnames(mod_full)) {
    mod <- mod_full[, !(colnames(mod_full) %in% reference_col)]
  } else {
    warning(paste("Reference timepoint column", reference_col, "not found. Using full model."))
    mod <- mod_full
  }
  
  # Run CovBat harmonization
  covbat_result <- covbat(
    dat = dat,
    bat = bat,
    mod = mod,
    percent.var = percent_var,
    std.var = standardize_var,
    mean.only = mean_only,
    eb = empirical_bayes,
    parametric = parametric,
    score.eb = FALSE,
    score.parametric = TRUE,
    resid = FALSE,
    verbose = verbose
  )
  
  # Extract harmonized data and transpose back
  harmonized_data <- t(covbat_result$dat.covbat)
  
  # Set column names with prefix
  colnames(harmonized_data) <- paste0("harmonisedCovbat_", feature_cols)
  
  # Combine with metadata
  # Only include metadata columns that exist in the data
  available_metadata <- metadata_cols[metadata_cols %in% colnames(df)]
  if (length(available_metadata) == 0) {
    warning("No metadata columns found. Returning harmonized data only.")
    final_df <- as.data.frame(harmonized_data)
  } else {
    final_df <- cbind(df[, available_metadata, drop = FALSE], harmonized_data)
  }
  
  # Return both the final dataframe and the full CovBat result object
  return(list(
    harmonized_df = final_df,
    covbat_result = covbat_result
  ))
}