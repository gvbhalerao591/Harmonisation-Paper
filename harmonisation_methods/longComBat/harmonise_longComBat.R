rm(list = ls())
harmonise_longComBat <- function(data, 
                                    batch_variable, 
                                    idp_names,
                                    age_var = "Final_Age",
                                    timepoint_var = "timepoint",
                                    subject_var = "subjectID",
                                    site_var = "Site",
                                    scanner_var = "Scanner",
                                    zscore_age = TRUE,
                                    save_output = FALSE,
                                    output_filename = NULL) {
  #' Harmonise Brain Volume Data Using longCombat
  #'
  #' @param data Data frame containing imaging data and covariates
  #' @param batch_variable Batch configuration: "site", "scanner", or "scanner_site"
  #'                       "site" = site as fixed effect
  #'                       "scanner" = scanner as fixed effect  
  #'                       "scanner_site" = scanner fixed + site random
  #' @param idp_names Character vector of imaging measure column names to harmonise
  #' @param age_var Name of age column (default: "Final_Age")
  #' @param timepoint_var Name of timepoint column (default: "timepoint")
  #' @param subject_var Name of subject ID column (default: "subjectID")
  #' @param site_var Name of site column (default: "Site")
  #' @param scanner_var Name of scanner column (default: "Scanner")
  #' @param zscore_age Logical, whether to z-score age (default: TRUE)
  #' @param save_output Logical, whether to save harmonised data to CSV (default: FALSE)
  #' @param output_filename Optional output filename (default: auto-generated)
  #'
  #' @return List containing:
  #'   - data_combat: harmonised data frame
  #'   - gammahat: location parameter estimates
  #'   - delta2hat: scale parameter estimates
  #'   - gammastarhat: empirical Bayes location estimates
  #'   - delta2starhat: empirical Bayes scale estimates
  #'
  #' @examples
  #' # Site as fixed effect
  #' result <- harmonise_brain_volumes(data, "site", c("vol1", "vol2"))
  #'
  #' # Scanner as fixed effect
  #' result <- harmonise_brain_volumes(data, "scanner", c("vol1", "vol2"))
  #'
  #' # Scanner fixed + site random
  #' result <- harmonise_brain_volumes(data, "scanner_site", c("vol1", "vol2"))
  #'
  #' # Custom variable names
  #' result <- harmonise_brain_volumes(data, "site", c("vol1"), 
  #'                                   age_var = "age_years",
  #'                                   site_var = "acquisition_site")
  
  # Load required libraries
  if (!require("longCombat", quietly = TRUE)) {
    stop("Package 'longCombat' is required. Install with: install.packages('longCombat')")
  }
  if (!require("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required. Install with: install.packages('lme4')")
  }
  
  # Validate inputs
  validate_inputs(data, batch_variable, idp_names, age_var, timepoint_var, 
                  subject_var, site_var, scanner_var)
  
  # Make a copy of data to avoid modifying original
  data_work <- data
  
  # Z-score age if requested
  if (zscore_age && age_var %in% names(data_work)) {
    if (is.numeric(data_work[[age_var]])) {
      data_work[[age_var]] <- scale(data_work[[age_var]])[, 1]
      message(paste0("Z-scored age variable: ", age_var))
    } else {
      warning(paste0("Age variable '", age_var, "' is not numeric. Skipping z-scoring."))
    }
  }
  
  # Select required columns
  required_cols <- unique(c(subject_var, timepoint_var, age_var, 
                            scanner_var, site_var, idp_names))
  required_cols <- intersect(required_cols, names(data_work))
  data_work <- data_work[, required_cols, drop = FALSE]
  
  # Convert variables to factors
  data_work[[site_var]] <- as.factor(data_work[[site_var]])
  data_work[[scanner_var]] <- as.factor(data_work[[scanner_var]])
  data_work[[subject_var]] <- as.factor(data_work[[subject_var]])
  data_work[[timepoint_var]] <- as.factor(data_work[[timepoint_var]])
  
  message(paste0("Found ", length(idp_names), " imaging measures to harmonise"))
  message(paste0("Found ", nrow(data_work), " observations"))
  message(paste0("Found ", nlevels(data_work[[subject_var]]), " subjects"))
  
  # Build formula and determine batch configuration
  formula_str <- paste(age_var, timepoint_var, sep = " + ")
  
  # Normalize batch_variable input
  batch_variable <- tolower(gsub("_", "", batch_variable))
  
  if (batch_variable == "site") {
    # Site as fixed effect
    batchvar_use <- site_var
    ranef_use <- paste0("(1|", subject_var, ")")
    outname <- "Site"
    message("Configuration: Site as FIXED batch effect, subject as random effect")
    
  } else if (batch_variable == "scanner") {
    # Scanner as fixed effect
    batchvar_use <- scanner_var
    ranef_use <- paste0("(1|", subject_var, ")")
    outname <- "Scanner"
    message("Configuration: Scanner as FIXED batch effect, subject as random effect")
    
  } else if (batch_variable == "scannersite") {
    # Scanner fixed + site random
    batchvar_use <- scanner_var
    ranef_use <- paste0("(1|", subject_var, ") + (1|", site_var, ")")
    outname <- "ScannerF_SiteR"
    message("Configuration: Scanner as FIXED batch effect, site as RANDOM batch effect, subject as random effect")
    
  } else {
    stop("batch_variable must be one of: 'site', 'scanner', or 'scanner_site'")
  }
  
  # Run longCombat
  message("\nRunning longCombat harmonisation...")
  result <- longCombat::longCombat(
    idvar = subject_var,
    timevar = timepoint_var,
    batchvar = batchvar_use,
    features = idp_names,
    formula = formula_str,
    ranef = ranef_use,
    data = data_work,
    verbose = TRUE
  )
  
  message("\nHarmonisation complete!")
  
  # Rename harmonised columns from ".combat" to "harmonisedLC_" prefix
  combat_cols <- paste0(idp_names, ".combat")
  new_cols <- paste0("harmonisedLC_", idp_names)
  
  # Find which columns exist (in case some failed)
  existing_combat_cols <- intersect(combat_cols, names(result$data_combat))
  
  if (length(existing_combat_cols) > 0) {
    # Get indices of combat columns
    combat_idx <- match(existing_combat_cols, names(result$data_combat))
    
    # Replace with new names
    for (i in seq_along(combat_idx)) {
      old_name <- existing_combat_cols[i]
      new_name <- gsub("\\.combat$", "", old_name)  # Remove .combat
      new_name <- paste0("harmonisedLC_", new_name)  # Add prefix
      names(result$data_combat)[combat_idx[i]] <- new_name
    }
    
    message(paste0("Renamed ", length(existing_combat_cols), 
                   " harmonised columns with 'harmonisedLC_' prefix"))
  }
  
  # Save output if requested
  if (save_output) {
    if (is.null(output_filename)) {
      output_filename <- paste0("harmonised_data_", outname, "_", 
                                format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    }
    write.csv(result$data_combat, output_filename, row.names = FALSE)
    message(paste0("Saved harmonised data to: ", output_filename))
  }
  
  return(result)
}


validate_inputs <- function(data, batch_variable, idp_names, age_var, 
                            timepoint_var, subject_var, site_var, scanner_var) {
  #' Validate function inputs
  
  # Check data
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  
  if (nrow(data) == 0) {
    stop("'data' is empty")
  }
  
  # Check batch_variable
  batch_variable <- tolower(gsub("_", "", batch_variable))
  valid_batch <- c("site", "scanner", "scannersite")
  if (!batch_variable %in% valid_batch) {
    stop("batch_variable must be one of: 'site', 'scanner', or 'scanner_site'")
  }
  
  # Check idp_names
  if (!is.character(idp_names) || length(idp_names) == 0) {
    stop("'idp_names' must be a non-empty character vector")
  }
  
  missing_idps <- setdiff(idp_names, names(data))
  if (length(missing_idps) > 0) {
    stop(paste0("The following IDPs are not in data: ", 
                paste(missing_idps, collapse = ", ")))
  }
  
  # Check required columns
  required_vars <- c(age_var, timepoint_var, subject_var)
  
  # Add batch-specific required columns
  if (batch_variable %in% c("site", "scannersite")) {
    required_vars <- c(required_vars, site_var)
  }
  if (batch_variable %in% c("scanner", "scannersite")) {
    required_vars <- c(required_vars, scanner_var)
  }
  
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(paste0("The following required variables are not in data: ", 
                paste(missing_vars, collapse = ", ")))
  }
  
  # Check for missing data in key variables
  check_vars <- c(required_vars, idp_names)
  missing_data <- sapply(data[, check_vars, drop = FALSE], function(x) sum(is.na(x)))
  if (any(missing_data > 0)) {
    vars_with_missing <- names(missing_data)[missing_data > 0]
    warning(paste0("Missing data detected in: ", 
                   paste(vars_with_missing, collapse = ", "), 
                   "\nlongCombat will fail if missing data exists. Consider imputation or removal."))
  }
  
  invisible(TRUE)
}


# Example usage demonstration
if (FALSE) {
  # Example 1: Site as fixed effect
  result1 <- harmonise_brain_volumes(
    data = my_data,
    batch_variable = "site",
    idp_names = c("T1_SIENAX_GM_norm_vol", "T1_SIENAX_WM_norm_vol")
  )
  
  # Example 2: Scanner as fixed effect with custom variable names
  result2 <- harmonise_brain_volumes(
    data = my_data,
    batch_variable = "scanner",
    idp_names = c("volume1", "volume2"),
    age_var = "age_years",
    scanner_var = "mri_scanner"
  )
  
  # Example 3: Scanner fixed + site random, save output
  result3 <- harmonise_brain_volumes(
    data = my_data,
    batch_variable = "scanner_site",
    idp_names = c("volume1", "volume2"),
    save_output = TRUE,
    output_filename = "harmonised_volumes.csv"
  )
  
  # Access harmonised data
  harmonised_data <- result3$data_combat
}