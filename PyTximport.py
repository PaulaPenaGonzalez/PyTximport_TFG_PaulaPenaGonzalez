def PytximportP (files,
                type,
                TxIn=True,
                TxOut=False,
                countsFromAbundance = None,
                tx2gene = None ,
                varReduce=False,
                dropInfReps=False,
                infRepStat=None,
                ignoreTxVersion=False,
                ignoreAfterBar=False,
                geneIdCol=None,
                txIdCol=None,
                abundanceCol=None,
                countsCol=None,
                lengthCol=None,
                importer=None,
                existenceOptional=False,
                sparse=False,
                sparseThreshold=1,
                readLength=75):

  infRepImporter = None,
  import sys
  import os
  import gzip
  import re
  import json
  import pandas as pd
  import numpy as np
  from scipy.sparse import csr_matrix

#Configuración del argumento type:
  if type.lower() not in ["kallisto", "salmon", "sailfish"]:
     raise ValueError("Packages other than Kallisto, Salmon, and Sailfish are not yet supported in PyTximport.")

#Configuración del argumento CountsFromAbundance
  if countsFromAbundance == "dtuScaledTPM" and not TxOut: #Si el argumento CountsFromAbundance es dtuScaledTMP, la normalización de los TMP se debe hacer a nivel de tránscritos y no a nivel de gen, se imprime un error
    raise ValueError("dtuScaledTMP can't be performed at gene-level: TxOut should be True")

#Configuración del argumento tx2gene (se asegura de que introduces un DataFrame o archivo .csv).
  if tx2gene is None:
    raise ValueError("tx2gene argument should be given as a two column DF or .csv file.")
  if not isinstance(tx2gene, pd.DataFrame):  # Si tx2gene no es un DataFrame:
        print("tx2gene input is not a DataFrame. Trying to open a .csv file...")
        if isinstance(tx2gene, str) and tx2gene.lower().endswith('.csv'):  #Verifica si es un string y si termina en .csv
            try:
                tx2gene = load_tx2gene(tx2gene)
            except Exception as e:
                raise ValueError(f"Failed to read the CSV file: {e}")
        else:
            raise ValueError("tx2gene must be a DataFrame or a path to a .csv file.")

  if tx2gene.shape[1] != 2:  #Verifica que el DataFrame tenga dos columnas:
    raise ValueError("tx2gene must be a 2 column DataFrame: first column for transcripts IDs and second one for genes IDs")


#Configuración del argumento existenceOptional que se asegura de que los archivos existen antes de intentar importarlos.
  if existenceOptional: #Si existenceOptional es True se va a comprobar la existencia de los archivos.
    for archivo in files:
      if os.path.exists(archivo):
        print(f"The file {archivo} exists.")
      else:
        print(f"The file {archivo} doesn't exist.")
        raise ValueError("Review the existence of the files listed in the list for proper importation.")

#Asegura que hay archivos en la lista files
  if not len(files)>0:
    raise ValueError("There are no files in the list")

#Configuración de los argumentos TxIn y TxOut: no es posible proporcionar una salida a nivel de transcripción si la entrada no lo es.
  if TxOut and not TxIn:
    raise ValueError("TxOut is only an option when transcript-level data is read in (TxIn=True)")

#Configuración para Kallisto
  if "abundance.h5" == os.path.basename(files[0]): #Si el archivo de abundancia está en formato .h5, la variable kallisto_h5 se guarda como True
    (kallisto_h5) = True
  else:
    (kallisto_h5) = False

  if type.lower() == "kallisto" and not kallisto_h5:
    raise ValueError("Importing 'abundance.h5' is faster than 'abundance.tsv'. Importing 'abundance.tsv' is not supported yet.") #Si no hay archivo de abundancia en .h5, se imprime un mensaje de error.

#The tximport arguments varReduce and dropInfReps can be used to summarize the inferential replicates into a single variance per transcript and per sample, or to not import inferential replicates, respectively.

#Configuración del argumento infRepStat
  if infRepStat is not None: #Si se ha introducido un valor para el argumento infRepStat se necesitan las Inf Reps y
    if dropInfReps:
      raise ValueError("You can't drop InfReps because infRepStat requires them")
    if type == "alevin":
      raise ValueError("infRepStat is presently incompatible with alevin input")
    if sparse:
      raise ValueError("infRepStat is not currently compatible with sparse output")
    if not callable(infRepStat):
      if infRepStat == "default":
        def infRepStat_function(matriz):
          medianas_por_fila = np.median(matriz, axis=1)
          return medianas_por_fila
        infRepStat = infRepStat_function
      else:
        raise TypeError("infRepStat must be a function if it is not 'default'.")


  #Configuración de readrStatus
  if not kallisto_h5 and importer == None:
    readrStatus = True #En R se establece como true cuando se usa el paquete readr para leer los archivos.

#Configuración columnas para Salmon y Sailfish
  typesS=["salmon","sailfish"]

  if type.lower() in typesS:
    txIdCol = "Name"
    abundanceCol = "TPM"
    countsCol = "NumReads"
    lengthCol = "EffectiveLength"
    if importer == None:
      def importer(x): #Definimos una función para importar los datos, donde X es cada archivo en files
        #Los tipos de objetos que esperamos encontrar en cada columna
        column_types = {"col1": str, "col2": int, "col3": float, "col4": float, "col5": float}
        data = pd.read_csv(x, sep="\t", dtype=column_types) #Leer los archivos con los tipos de columna especificados

      #Si se quieren leer o no las réplicas inferenciales de los datos de las muestras
      if dropInfReps: #Si True: no se quieren leer: infRepImporter vacío
        infRepImporter = None #NO creamos función para importar réplicas inferenciales
      else:
        infRepImporter = readInfRepFish #Igualamos a función predefinida para importar las réplicas inferenciales


#Configuración para Kallisto:

  if type.lower() == "kallisto":
    txIdCol = "target_id"
    abundanceCol = "tpm"
    countsCol = "est_counts"
    lengthCol = "eff_length"
    if kallisto_h5:
      importer = read_kallisto_h5
    elif not kallisto_h5 and importer == None:
      def importer(x): #Definimos una función para importar los datos, donde X es cada archivo en files
          #Los tipos de objetos que esperamos encontrar en cada columna
          column_types = {"col1": str, "col2": int, "col3": float, "col4": float, "col5": float}
          data = pd.read_csv(x, sep="\t", dtype=column_types)

    #Si se quieren leer o no las réplicas inferenciales de los datos de las muestras
    if dropInfReps: #Si True: no se quieren leer: infRepImporter vacío
      infRepImporter = None #No se crea función para importar réplicas inferenciales
    else:
      infRepImporter = readInfRepKallisto #Igualamos a función predefinida para importar las réplicas inferenciales

#InfRepType inicialización:
  infRepType = "none"
  #Si se quieren importar las réplicas inferenciales:
  if type.lower() in ["kallisto", "salmon", "sailfish"] and not dropInfReps:
    if varReduce and TxOut: #Solo varianza de las réplicas.
      infRepType = "var"
    else:
      infRepType = "full"

  #Prueba de la importacion de réplicas inferenciales:
  repInfo=None
  if infRepType != "none":
    if type.lower() in typesS:
      repInfo = infRepImporter(os.path.dirname(files[0]))
    else:
      repInfo=None
      print("Only Salmon and Sailfish types are available")

    if repInfo == None:
      infRepType = "none"

  # Condiciones para importación de los datos dispersos: de matriz dispersa a matriz densa
  if sparse:
    SvalidCFA = ["no", "scaledTPM"]
    if not infRepType == "none" and not TxOut and countsFromAbundance not in SvalidCFA:
      raise ValueError("Importing sparsely, only counts and abundace would be returned. Support is restricted to scenarios where txOut is TRUE, CFA is either 'no' or 'scaledTPM', and there are no inferential replicates ")

  ## ------- Bucle principal -------
  for i in range(len(files)):
    print(f"Iteration over file {i}:")
    raw = importer(files[i]) #Función que devuelve un df de Pandas
    if infRepType != "none":
      if type.lower() == "piscem":
        repInfo = infRepImporter(files[i])

      #if type.lower() in ["salmon","sailfish", ]:
      else:
        repInfo = infRepImporter(os.path.dirname(files[i]))

    else:
      repInfo = None
      print("Inferential replicates will not be imported.")

    if repInfo == None:
      infRepType = "none"
      print("RepInfo is none. Couldn't import InfReps from incompatible types.")

    #Comprobar que están todas las columnas en el DF
    if not all(col in raw.columns for col in [abundanceCol, countsCol, lengthCol]):
      raise ValueError("Not all required columns are present in the DataFrame.")

    #Comprobación cruzada de los valores de las columnas
    if i==0:
      txId=raw[txIdCol]
    else:
      if not raw[txIdCol].equals(txId):
        raise ValueError("Los valores en la columna txId no son consistentes entre las iteraciones. (Prueba)")
      else:
        print("Los valores de las columnas son consistentes entre muestras. (Prueba)")

    #Importar matrices densas desde matrices dispersas
    if not sparse:
      if i == 0:
        #Como es la primera iteración creamos matrices vacías con tantas filas como tiene raw y tantas columnas como archivos.
        #Guardamos los nombres de filas y columnas en variables (no podemos guardarlas en el array de numpy)
        matrix_row_names = raw[txIdCol].values #Nombres de las filas como los de la columna txId de raw
        matrix_col_names = files #Nombres de los archivos como nombres de las columnas
        #Asignamos la matriz vacía a las variables para que tengan todas el mismo tamaño.
        abundance_mat_tx = np.empty((raw.shape[0], len(files)))
        counts_mat_tx = np.empty((raw.shape[0], len(files)))
        length_mat_tx = np.empty((raw.shape[0], len(files)))
        if infRepType == "var":
          var_mat_tx = np.empty((raw.shape[0], len(files)))
        if infRepType == "full":
          inf_rep_dict_tx = {} #Se asigna un diccionario vacío

      #Comenzamos a llenar las matrices vacías con los valores
      abundance_mat_tx[:, i] = raw[abundanceCol]
      counts_mat_tx[:, i] = raw[countsCol]
      length_mat_tx[:, i] = raw[lengthCol]
      if infRepType == "var":
        var_mat_tx[:, i] = repInfo["vars"]
      if infRepType == "full":
        #Se crea un matriz de una sola fila (1D) con todos los valores de "reps" y se guarda en su clave i del diccionario
        inf_rep_dict_tx[i] = repInfo["reps"] #Se guarda la matriz de reps en su clave "i" del diccionario

      #Configuración con infRepStat predefinido como argumento de la función
      if infRepStat != None:
        counts_mat_tx[:, i] = infRepStat(repInfo["reps"]) #Donde infRepStat es una función predefinida que vuelve a calcular conteos y abundancias a partir de las réplicas inferenciales, sacando la mediana.
        tpm = counts_mat_tx[:, i] / length_mat_tx[:, i]
        abundance_mat_tx[:,i] = tpm * 1e6 / sum(tpm) #Esto normaliza las abundancias de transcripción para que la suma total sea igual a 1 millón.

    else: #Si queremos quitar la dispersión de las matrices: sparse=True.
      sparse_idx = raw.index[raw[countsCol] >= sparseThreshold]

      if i == 0:
        txId = raw[txIdCol]
        NumNonZeros=[]
        counts_Array_I=np.empty((1, len(sparse_idx))).flatten()
        counts_Array_X=np.empty((1, len(sparse_idx))).flatten()

      if not raw[txIdCol].equals(txId):
        raise ValueError("Los valores en la columna txId no son consistentes entre las iteraciones.")
      else:
        print("Los valores de las columnas son consistentes entre muestras")

      NumNonZeros.append(len(sparse_idx))

      #Concatenamos en un array las abundancias, el primer elemento del array es un 0 por tanto la indexación ahora comienza en 1
      if i ==0:
        counts_Array_I[:]= np.array(sparse_idx)
        counts_Array_X[:]=raw.loc[counts_Array_I, countsCol].to_numpy()
      else:
        counts_Array_X = np.concatenate((counts_Array_X, raw.loc[sparse_idx, countsCol].to_numpy())) #Selecciona los elementos de la columna countsCol en las filas correspondientes a los índices almacenados en sparse_idx.
        counts_Array_I = (np.concatenate((counts_Array_I, np.array(sparse_idx)))).astype(int)

      #Si queremos hacer normalización por TPMs, sacamos las abundancias:
      if countsFromAbundance.lower() == "scaledtpm":
        if i==0:
          abundance_Array_X = np.empty((1,len(sparse_idx) )).flatten()
          abundance_Array_X[:]=raw.loc[counts_Array_I, abundanceCol].to_numpy()
          print("Array de abundancias de las cuentas no dispersas:",abundance_Array_X,
              np.shape(abundance_Array_X))
        else:
          #Concatenamos en un array las abundancias, el primer elemento del array es un 0 por tanto la indexación ahora comienza en 1
          abundance_Array_X = np.concatenate((abundance_Array_X, raw.loc[sparse_idx, abundanceCol].to_numpy()))
          print("Array de abundancias de las cuentas no dispersas:",abundance_Array_X,
              np.shape(abundance_Array_X))

  #----Fin Bucle Principal----

  if sparse:
    #Creamos una matriz dispersa con tantas filas como len(txId) y tantas columnas como muestras
    counts_mat_tx = csr_matrix((counts_Array_X, (counts_Array_I, np.repeat(np.arange(len(NumNonZeros)), NumNonZeros))), shape=(len(txId), len(files)))
    print(counts_mat_tx,
          np.shape(counts_mat_tx))
    if countsFromAbundance.lower() == "scaledtpm":
      abundance_mat_tx = csr_matrix((abundance_Array_X,(counts_Array_I, np.repeat(np.arange(len(NumNonZeros)), NumNonZeros))), shape=(len(txId), len(files)))

    else:
      abundance_mat_tx = None

    length_mat_tx = None

  # propagate names to inferential replicate list
  if infRepType == "full":
    if len(inf_rep_dict_tx) == len(files):
      infRepDictTx = {}
      for nombre, valor in zip(files, inf_rep_dict_tx.values()):
        infRepDictTx[nombre] = valor
      print(infRepDictTx)
    else:
      raise ValueError("Note: not all samples contain inferential replicates. Tximport can only import data when either all or no samplescontain inferential replicates. Instead first subset to theset of samples that all contain inferential replicates.")
  print("")

  #Creamos un bucle para sacar los nombres de las muestras y ponerlos como nombres de las columnas del DF
  folder_names = []

  # Obtener el nombre de la carpeta de cada archivo y agregarlo a la lista
  for file in files:
    folder_name = os.path.basename(os.path.dirname(file))
    folder_names.append(folder_name)

  abundance_df = pd.DataFrame(abundance_mat_tx, index=txId, columns=folder_names)
  counts_df = pd.DataFrame(counts_mat_tx, index=txId, columns=folder_names)
  length_df = pd.DataFrame(length_mat_tx, index=txId, columns=folder_names)

  if infRepType == "none":
    txi = {"abundance":abundance_df,
                "counts":counts_df,
                "length":length_df,
                "countsFromAbundance":countsFromAbundance}

  elif infRepType == "var":
    var_df = pd.DataFrame(var_mat_tx, index=txId, columns=folder_names)
    #Si solo mantenemos la varianza de las réplicas inferenciales
    txi = {"abundance":abundance_df,
                "counts":counts_df,
                "variance":var_df,
                "length":length_df,
                "countsFromAbundance":countsFromAbundance}

  elif infRepType == "full":
    #Si mantenemos todas las muestras de las réplicas inferenciales
    txi = {"abundance":abundance_df,
           "counts":counts_df,
           "infReps":inf_rep_dict_tx,
           "length":length_df,
           "countsFromAbundance":countsFromAbundance}


  #Dos tipos de outputs basados en la configuración de TxOut:

  #1) Si el usuario quiere los datos a nivel de transcripción:
  if TxOut:
    if countsFromAbundance.lower() != "no":
      #Para dtuScaledTPM, actuar como si fuese lengthScaledTPM pero con una matriz de longitudes alterada.
      length4CFA = txi["length"] # intermediate version of the length matrix
      if countsFromAbundance.lower() == "dtuscaledtpm":
        length4CFA = medianLengthOverIsoform(length4CFA, tx2gene)
        countsFromAbundance = "lengthScaledTPM"
      #Función para calcular cualquiera de los 3 métodos de countsFromAbundance:
      txi["counts"] = makeCountsFromAbundance(txi["counts"],
                                              txi["abundance"],
                                              length4CFA,
                                              countsFromAbundance)
    # Para cambiar los nombres de las filas hbaría que transformar las matrices a DF de pandas:
    # Código de R comentado para traducir si hace falta
    #if ignoreAfterBar:
      #for (matNm in c("counts","abundance","length")):
        #rowNms <- rownames(txi[[matNm]])
        #rownames(txi[[matNm]]) <- sub("\\|.*", "", rowNms)
    return(txi)

  #2)Por el contrario, si se quiere la salida a nivel de gen:
  else:
    txi["countsFromAbundance"] = None
    txiGene = summarizeToGene(txi, tx2gene, varReduce, ignoreTxVersion, ignoreAfterBar, countsFromAbundance)

    return(txiGene)