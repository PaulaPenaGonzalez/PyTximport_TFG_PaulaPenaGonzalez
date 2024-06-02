import numpy as np
from scipy.sparse import csr_matrix, isspmatrix_csr

def makeCountsFromAbundance(countsMatrix, abundanceMatrix, lengthMat, countsFromAbundance):
  PossibleOptions = ["scaledtpm","lengthscaledtpm"]

  if countsFromAbundance.lower() in PossibleOptions:
    # Comprobar que la matriz de cuentas es de tipo CSR
    sparse = isspmatrix_csr(countsMatrix) #Innecesario si vamos a usar la misma función sea matriz densa o dispersa

    # En Python usamos la misma función sum(axis=0) para hacer la suma de las columnas tanto en matrices densas como en dispersas
    def sumColFun(countsMatrix):
      #countsSum = countsMatrix.sum(axis=0)
      countsSum = np.nansum(countsMatrix, axis=0)
      return countsSum

    # Calcular el sumatorio de las columnas de la matriz de cuentas
    countSum = sumColFun(countsMatrix)

    # En función del tipo de normalización:
    if countsFromAbundance.lower() == "lengthscaledtpm":
      # Calcular nuevas cuentas:
      newCounts = (abundanceMatrix * (np.mean(lengthMat, axis=1))[:, np.newaxis])

    elif countsFromAbundance.lower() == "scaledtpm":
      newCounts = abundanceMatrix

    else:
      raise ValueError("expecting 'lengthScaledTPM' or 'scaledTPM'")

    newSum = sumColFun(newCounts)

    if sparse:
      # Hay que probar si funciona con matrices dispersas de csr
      countsMat = np.transpose(np.transpose(newCounts)*((countSum / newSum)[:, np.newaxis]))
    else:
      countsMat=np.transpose(np.transpose(newCounts)*((countSum / newSum)[:, np.newaxis]))

  return countsMat