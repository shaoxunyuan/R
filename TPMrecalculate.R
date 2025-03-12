TPMrecalculate <- function(quant.sf.data) {
  # Check if the input data frame contains "sampleID" or "BioSample" columns
  # If neither is found, stop execution and return an error message
  if (!("sampleID" %in% names(quant.sf.data) || "BioSample" %in% names(quant.sf.data))) {
    stop("Neither 'sampleID' nor 'BioSample' found in quant.sf.data.")
  }
  # Determine which column to use as the sample identifier column
  # If the data frame contains "sampleID" column, use "sampleID" as the sample identifier column
  # Otherwise, use "BioSample" as the sample identifier column
  sample_col <- ifelse("sampleID" %in% names(quant.sf.data), "sampleID", "BioSample")
  # Print the current sample identifier column being used
  print(paste("Using column:", sample_col))
  # Print message indicating that we are checking the filtered columns: sample identifier column, number of reads column, and effective length column
  print("Checking filtered columns: sampleID/BioSample, NumReads, EffectiveLength")
  # Initialize an empty list to store the recalculated TPM values for each sample
  TPM.adj.list <- list()
  # Loop through each unique sample in the data frame
  for (onesample in unique(quant.sf.data[[sample_col]])) {
    # Extract the data for the current sample
    quant.sf.data.onesample <- quant.sf.data[quant.sf.data[[sample_col]] == onesample, ]
    # Recalculate the TPM values for the current sample
    # TPM calculation formula: (each transcript's number of reads / that transcript's effective length) / (sum of all transcripts' number of reads / effective length) * 1000000
    quant.sf.data.onesample$TPMadj <- (quant.sf.data.onesample$NumReads / quant.sf.data.onesample$EffectiveLength) /
                                      (sum(quant.sf.data.onesample$NumReads / quant.sf.data.onesample$EffectiveLength)) * 1000000
    # Store the recalculated data for the current sample in the list
    TPM.adj.list[[onesample]] <- quant.sf.data.onesample
  }
  # Combine all data frames in the list into one large data frame by rows
  TPM.adj.df <- do.call(rbind, TPM.adj.list)
  # Convert the combined data frame to data frame type (this step is redundant and can be removed)
  # TPM.adj.df <- as.data.frame(TPM.adj.df)
  # Reset the row names of the data frame
  rownames(TPM.adj.df) <- NULL
  # Return the data frame with recalculated TPM values
  return(TPM.adj.df)
}
