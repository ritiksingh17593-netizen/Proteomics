library(Spectra)
library(magrittr)

file <- "C:/Users/LENOVO/Desktop/Proteomics/data/raw/PXD062725/DIA_iMN_WT_2uL_r5.mzML.xml"

sps <- Spectra(file) 
table(msLevel(sps))


i1 <- which(msLevel(sps) == 1)[1]
i2 <- which(msLevel(sps) == 2)[1]
i1
i2

# Extract the peak matrix for the first MS1 spectrum before plotting it.
ms1_peaks <- peaksData(sps[i1])[[1]]

spectraData(
    sps[i1],
    c("msLevel", "rtime", "scanIndex")
)

spectraData(
    sps[i2],
    c("msLevel", "rtime", "scanIndex")
)

plot(
    ms1_peaks[,1],
    ms1_peaks[,2],
    type = "h",
    xlab = "m/z",
    ylab = "Intensity",
    main = "First MS1 Spectrum"
)

plot(
    ms2_peaks[, 1],
    ms2_peaks[, 2],
    type = "h",
    xlab = "m/z",
    ylab = "Intensity",
    main = "First MS2 Spectrum"
)


p2 <- peaksData(sps[i2])[[1]]

plot(
    p2[, 1],
    p2[, 2],
    type = "h",
    xlab = "m/z",
    ylab = "Intensity",
    main = "First MS2 Spectrum"
)

