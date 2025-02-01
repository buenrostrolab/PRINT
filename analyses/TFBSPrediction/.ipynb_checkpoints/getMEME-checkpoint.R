# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# First load cisBP motifs in PWMatrix format
species <- "mouse"
cisBPMotifs <- readRDS(paste0("/data/PRINT/multiScaleFootprinting/data/shared/cisBP_", species, "_pwms_2021.rds"))

# Convert to MEME format
outputPath <- paste0("/data/PRINT/multiScaleFootprinting/data/shared/cisBP_", species, "_pwms_2021.meme")

# Write header lines
header <- c("MEME version 4\n", "ALPHABET= ACGT\n", "strands: + -",
            "Background letter frequencies", "A 0.25 C 0.25 G 0.25 T 0.25\n")
for(h in header){
  write(h, file = outputPath, append = TRUE)
}

# Write PPMs of each TF
for(TF in names(cisBPMotifs)){
  print(paste(TF, "Progress", match(TF, names(cisBPMotifs))))
  PPM <- (2 ^ cisBPMotifs[[TF]]@profileMatrix) / 4
  write(paste("MOTIF", TF), file = outputPath, append = TRUE)
  write(paste("letter-probability matrix: alength= 4 w=", dim(PPM)[2], " nsites= 200"), file = outputPath, append = TRUE)
  for(i in 1:dim(PPM)[2]){
    write(paste(sprintf("%.6f", round(PPM[, i], 6)), collapse = " "), file = outputPath, append = TRUE)
  }
  write("", file = outputPath, append = TRUE)
}
