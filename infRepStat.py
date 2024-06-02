def infRepStat(matriz):
  import numpy as np
  medianas_por_fila = np.median(matriz, axis=1)
  return medianas_por_fila
