- TFBSTrainingData.R: This script loads the pre-computed multi-scale footprints at specified scales, as well as the GRanges object of ChIP-seq TF binding events (ENCODE or Unibind). Generates paired footprint and binding labels for later model training. This is saved to files like TFBSDataUnibind.h5

- footprint_to_TF.ipynb: Trains the footprint-to-TF prediction model which uses multi-scale footprints as input and predicts TF binding. Uses TFBSDataUnibind.h5 as input. Also contains code for generating the training data for seq2PRINT TF binding score (***_pred_data.tsv)

- 
