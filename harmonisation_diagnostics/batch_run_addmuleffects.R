# -------------------------------
# Self-mapping long_combat runner
# -------------------------------

raw_file_path <- "/Users/psyc1586_admin/GVB_data/harmonisation_work/00_aa_skullImprovement_Oct2025/00_running_for_paper/run_eval_metrics_IQMscannerOnly//data/addmuleffects/Raw.csv"
harmonised_dir <- "/Users/psyc1586_admin/GVB_data/harmonisation_work/00_aa_skullImprovement_Oct2025/00_running_for_paper/run_eval_metrics_IQMscannerOnly/data/addmuleffects"

setwd(harmonised_dir)

# Define the mapping: file names -> batch variable & formula
mapping <- list(
  "ComBat_Scanner.csv"           = list(batch="Site",   formula="zscore_age + timepoint"),
  "ComBat_Site.csv"              = list(batch="Site",      formula="zscore_age + timepoint"),
  "CovBat_Scanner.csv"           = list(batch="Site",   formula="zscore_age + timepoint"),
  "CovBat_Site.csv"              = list(batch="Site",      formula="zscore_age + timepoint"),
  "HistMatch.csv"                = list(batch="Site",      formula="zscore_age + timepoint"),
  "LME_IQM_Site.csv"             = list(batch="Site",      formula="zscore_age + timepoint"),
  "LME_IQM_Scanner.csv"          = list(batch="Site",   formula="zscore_age + timepoint"),
  "LME_Scanner.csv"              = list(batch="Site",   formula="zscore_age + timepoint"),
  "LME_ScannerFSiteR.csv"        = list(batch="Site",   formula="zscore_age + timepoint"),
  "LME_Site.csv"                 = list(batch="Site",      formula="zscore_age + timepoint"),
  "LongComBat_Scanner.csv"       = list(batch="Site",   formula="zscore_age + timepoint"),
  "LongComBat_ScannerFSiteR.csv" = list(batch="Site",   formula="zscore_age + timepoint"),
  "LongComBat_Site.csv"          = list(batch="Site",      formula="zscore_age + timepoint"),
  "LR_Scanner.csv"               = list(batch="Site",   formula="zscore_age + timepoint"),
  "LR_Site.csv"                  = list(batch="Site",      formula="zscore_age + timepoint"),
  "SynthSR.csv"                  = list(batch="Site",      formula="zscore_age + timepoint")
)

# Get all CSV files in harmonised_dir
all_files <- list.files(harmonised_dir, pattern="\\.csv$", full.names=FALSE)

# Filter only files that are in the mapping
files_to_run <- intersect(all_files, names(mapping))

# Loop through each file safely
for (file_name in files_to_run) {
  harmonised_file_path <- file.path(harmonised_dir, file_name)
  batch_var <- mapping[[file_name]]$batch
  formula_str <- mapping[[file_name]]$formula
  outdir <- file.path(harmonised_dir, tools::file_path_sans_ext(file_name))
  
  cat("Processing file:", file_name, "\n")
  
  run_AddMulEffects(
    raw_file_path = raw_file_path,
    harmonised_file_path = harmonised_file_path,
    outdir = outdir,
    batch_var = batch_var,
    formula_str = formula_str
  )
  
  cat("Finished file:", file_name, "\n\n")
}

cat("All files processed successfully!\n")
