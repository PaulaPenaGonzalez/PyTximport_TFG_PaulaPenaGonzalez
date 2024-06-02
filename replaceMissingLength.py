def replaceMissingLength(lengthMat, aveLengthSampGene):
  import numpy as np
  nanRows = lengthMat.isna().any(axis=1)
  if nanRows.any():
    nanRows.index = range(len(nanRows))
    indices_nanRows = nanRows[nanRows].index

  for index in indices_nanRows:
    if (lengthMat.iloc[index].isnull()).all():
      lengthMat.iloc[index] = aveLengthSampGene.iloc[index]

    else:
      nan_positions = lengthMat.iloc[index].isna()
      notNan_values = lengthMat.iloc[index][lengthMat.iloc[index].notna()]
      log_values = notNan_values.apply(np.log)
      median_log_values = np.mean(log_values)
      exponencial=np.exp(median_log_values)
      lengthMat.iloc[index][nan_positions] = exponencial

  return(lengthMat)