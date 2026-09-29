library(Spectra)

# 1. Load just one of your files for plotting
file_path <- file.path("data", "raw", "PXD062725", "DIA_iMN_dNLS_2uL_r1.mzML.xml")
sps <- Spectra(file_path)

# 2. Find the index number of the first MS1 scan and the first MS2 scan
i1 <- which(msLevel(sps) == 1)[1]
i2 <- which(msLevel(sps) == 2)[1]

# 3. Extract the peak matrices (m/z and intensity values)
# peaksData returns a list of matrices; [[1]] pulls the matrix for the scan
ms1_peaks <- peaksData(sps[i1])[[1]]
ms2_peaks <- peaksData(sps[i2])[[1]]

# 4. Plot the MS1 spectrum
plot(
  x = ms1_peaks[, 1],       # The first column is the m/z (mass-to-charge) ratio
  y = ms1_peaks[, 2],       # The second column is the Intensity
  type = "h",               # "h" draws vertical lines (needles) for mass spec peaks
  xlab = "m/z", 
  ylab = "Intensity", 
  main = "First MS1 Spectrum"
)

# 5. Plot the MS2 spectrum
plot(
  x = ms2_peaks[, 1], 
  y = ms2_peaks[, 2], 
  type = "h", 
  xlab = "m/z", 
  ylab = "Intensity", 
  main = "First MS2 Spectrum (Fragmentation)"
)