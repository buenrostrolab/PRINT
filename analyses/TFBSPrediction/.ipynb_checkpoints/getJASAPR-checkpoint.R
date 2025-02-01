# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
library(TFBSTools)

####################
# Load JASPAR file #
####################

path <- "/n/holylfs05/LABS/buenrostro_lab/Users/yanhu/cell_dynamics/multiScaleFootprinting/data/shared/JASPAR_2022.jaspar"
lines <- readLines(path) 
PWMList <- list()
for(line in lines){
  if(substr(line, 1, 1) == ">"){
    TF <- strsplit(line, "\t")[[1]][2]
    counter <- 0
    PFM <- NULL
  }else{
    counter <- counter + 1
    line_split <- strsplit(line, "[ \\[]")[[1]]
    line_split <- line_split[line_split != ""]
    PFM <- rbind(PFM, as.integer(line_split[2:(length(line_split) - 1)]))
    PFM <- PFM + 1 # Add pseudo-count
    if(counter == 4){
      PWM <- log2(t(t(PFM) / colSums(PFM)) / 0.25)
      rownames(PWM) <- c("A", "C", "G", "T")
      PWM <- TFBSTools::PWMatrix(name = TF, profileMatrix = as.matrix(PWM))
      PWMList[[TF]] <- PWM
    }
  }
}

## Construction of PFM<atrixList from list of PFMatrix
PWMList <- do.call(PWMatrixList, PWMList)

saveRDS(PWMList, "../../data/shared/JASPAR_human_pwms_2022.rds")
