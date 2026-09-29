library(Spectra)
library(dplyr)
library(tibble)

setwd("C:/Users/LENOVO/Desktop/Proteomics")

# 1. Locate all mzML files using a relative path
raw_dir <- file.path("data", "raw", "PXD062725")
mzml_files <- list.files(raw_dir, pattern = "mzML", full.names = TRUE)
message(paste("Found", length(mzml_files), "mzML files."))

# 2. Iterate over each file and extract QC metrics
qc_summary_list <- lapply(mzml_files, function(f) {
  fname <- basename(f)
  message("Processing: ", fname)
  
  # Load the file lazily (metadata only, keeps memory low)
  sps <- Spectra(f)
  
  # Scan levels
  levels <- msLevel(sps)
  ms1_n <- sum(levels == 1)
  ms2_n <- sum(levels == 2)
  total_scans <- length(sps)
  
  # Retention time in minutes
  rts <- rtime(sps) / 60
  rt_min <- round(min(rts, na.rm = TRUE), 2)
  rt_max <- round(max(rts, na.rm = TRUE), 2)
  
  tibble(
    file_name = fname,
    total_spectra = total_scans,
    ms1_scans = ms1_n,
    ms2_scans = ms2_n,
    rt_start_min = rt_min,
    rt_end_min = rt_max
  )
})

# Combine into a single QC table
qc_table <- bind_rows(qc_summary_list)
print(qc_table)

# 3. Export QC metrics to your results folder
dir.create("results", showWarnings = FALSE)
write.csv(qc_table, file.path("results", "raw_files_qc_summary.csv"), row.names = FALSE)


BiocManager::install("Spectra")
