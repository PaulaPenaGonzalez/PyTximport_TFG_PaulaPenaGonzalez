def summarizeToGene(txi, tx2gene, countsFromAbundance, varReduce=False, ignoreTxVersion=False, ignoreAfterBar=False):

  import re
  txId = raw[txIdCol]

  if countsFromAbundance.lower() not in ["no","scaledtpm","lengthscaledtpm"]:
    raise ValueError("countsFromAbundance must be set as 'no','scaledTPM'or'lengthScaledTPM' for gene summarization")
  if txi["countsFromAbundance"] != None:
    if countsFromAbundance.lower() == "no" and txi["countsFromAbundance"] != "no":
      warning_message = ("Los recuentos originales ya no son accesibles, ya que countsFromAbundance tiene el valor:", txi["countsFromAbundance"],
                         "Para utilizar countsFromAbundance='no', vuelve a ejecutar objectmport() con esta configuración. ",
                         "Anulando 'countsFromAbundance' para establecerlo en:", txi["countsFromAbundance"])
      print(warning_message)
      countsFromAbundance = txi["countsFromAbundance"]

  # unpack matrices from list for cleaner code
  abundanceMatTx = txi["abundance"]
  countsMatTx = txi["counts"]
  lengthMatTx = txi["length"]

  #listIds=(abundanceMatTx.index).tolist()
  IdsCounts=(countsMatTx.index).tolist()
  IdsLength=(lengthMatTx.index).tolist()

  #dfIds = pd.DataFrame(listIds, columns=['Ids de los tránscritos'])

  #if not ((countsMatTx.index).tolist()).equals(listIds) and not ((lengthMatTx.index).tolist()).equals(listIds):
  if not txId.isin(IdsCounts).all() and not txId.isin(IdsLength).all() :
        raise ValueError("Los valores de los Ids no son consistentes entre los df.")

  if not tx2gene.empty:
    if ignoreTxVersion:
      #Eliminamos todas las versiones detrás del punto para sacar las isoformas de los tránscritos.
      txId = txId.apply(lambda x: re.sub("\\..*", "", str(x)))

      tx2gene = tx2gene.set_axis(["tx","gene"], axis='columns')
      #Creamos una serie de pandas sin isoformas de los tránscritos
      tx2gene_tx = tx2gene["tx"].apply(lambda x: re.sub("\\..*", "", str(x)))

      # Crear una copia del DataFrame original y asignar los valores sin isoformas de la columna 'tx'
      tx2gene_modified = tx2gene.copy()
      tx2gene_modified['tx'] = tx2gene_tx

      # Cambiar el nombre de las columnas del DF de las anotaciones
      tx2gene = cleanTx2gene(tx2gene_modified)

    elif ignoreAfterBar:
      #Eliminamos las versiones detrás de la barra para sacar las isoformas.
      txId = txId.apply(lambda x: re.sub("\\|.*", "", str(x)))

    txId_df = pd.DataFrame({"trásncritosId":txId})

    # Cambiar el nombre de las columnas del DF de las anotaciones
    tx2gene = tx2gene.set_axis(["tx","gene"], axis='columns')
    tx2gene = cleanTx2gene(tx2gene)

    if not txId.isin(tx2gene["tx"]).all():
      #tx2gene_tx.isin(txId).all():
      print("Not all values of raw[txIdCol] column are present in 'tx' column of tx2gene.",
            "There will be transcrips with NaN gene anotation.")
    else:
      print("All values os Ids are present in tx2gene")


    #Calculamos la intersección entre txId y tx2gene para sacar los tráscritos anotados presentes
    #interseccion = pd.Series(list(set(tx2gene_cleaned['tx']) and set(txId)))
    interseccion = set(tx2gene['tx']).intersection(set(txId))
    interseccion = pd.Series(list(interseccion))

    #Creamos un df con la intersección:
    interseccion_df = pd.DataFrame({'elemento': interseccion})

    #Creamos el df con las filas con intersección: conservamos solo los genes que corresponden con interseccion_df
    merged_df = pd.merge(txId_df, tx2gene, left_on= 'trásncritosId', right_on='tx', how='left')

    # Extraer las filas con valores NaN en alguna columna
    filas_con_nan = merged_df[merged_df.isna().any(axis=1)]

    #Se limpia el df de tránscritos con anotaciones de gen == NaN (puede ser redundante)
    df_cleaned = merged_df.dropna()

    indices_lista = df_cleaned.index.tolist()

    # Agrupar los índices por el valor de la columna 'gene'
    grupos = df_cleaned.groupby('gene').groups
    longitud=len(grupos.keys())

    # Iterar sobre cada lista en el diccionario y convertir cada elemento a entero
    for clave, lista in grupos.items():
      grupos[clave] = [int(elemento) for elemento in lista]

    # Como sobran pares clave-valor (están vacías de índices), crear una lista de claves a eliminar
    claves_a_eliminar = []

    # Iterar sobre el diccionario
    for clave, valor in grupos.items():
      if not valor:
          claves_a_eliminar.append(clave)

    longitud_eliminar=len(claves_a_eliminar)

      # Eliminar las claves del diccionario
    for clave in claves_a_eliminar:
      del grupos[clave]

    longitud_nueva_grupos = len(grupos.keys())

    num_total_elementos = sum(len(lst) for lst in grupos.values())

    if num_total_elementos == len(indices_lista):
      print("La extracción de los índices de cada agrupación por gen se ha hecho correctamente.")

    indices_agrupados = list(grupos.values())

    # Calcular la suma por filas de las abundancias utilizando gene_indices
    abundanceMat_g = pd.DataFrame([abundanceMatTx.iloc[indices].sum() for indices in indices_agrupados], index=grupos.keys())
    abundanceMat_g_np = abundanceMat_g.values

    # Calcular la suma por filas de las cuentas utilizando gene_indices
    countsMat_g = pd.DataFrame([countsMatTx.iloc[indices].sum() for indices in indices_agrupados], index=grupos.keys())
    countsMat_g_np = countsMat_g.values


    # Inicializar el diccionario infReps_gDict
    infReps_gDict = {}

    if "infReps" in txi:
        infReps_dict = txi["infReps"]
        print("Summarizing inferential replicates")

        for key, values in infReps_dict.items():
            infReps_i = pd.DataFrame(infReps_dict[key])

            # Calcular la suma agrupada por genes
            infReps_g = pd.DataFrame([infReps_i.iloc[indices].sum() for indices in indices_agrupados], index=grupos.keys())

            # Almacenar el resultado en infReps_gDict
            infReps_gDict[key] = infReps_g


    # the next lines calculate a weighted average of transcript length,
    # weighting by transcript abundance.
    # this can be used as an offset / normalization factor which removes length bias
    # for the differential analysis of estimated counts summarized at the gene level.
    wheighted_length_tx = abundanceMatTx * lengthMatTx
    wheighted_length_g = pd.DataFrame([wheighted_length_tx.iloc[indices].sum() for indices in indices_agrupados], index=grupos.keys())
    lengthMat = wheighted_length_g / abundanceMat_g

    # pre-calculate a simple average transcript length
    # for the case the abundances are all zero for all samples.
    # first, average the tx lengths over samples
    aveLengthSamp = lengthMatTx.mean(axis=1)
    # then simple average of lengths within genes (not weighted by abundance)
    aveLengthSamp_g = pd.DataFrame([aveLengthSamp.iloc[indices].mean(axis=0) for indices in indices_agrupados], index=grupos.keys())


    # check for NaN and if possible replace these values with geometric mean of other samples.
    # (the geometic mean here implies an offset of 0 on the log scale)
    # NaN come from samples which have abundance of 0 for all isoforms of a gene, and
    # so we cannot calculate the weighted average. our best guess is to use the average
    # transcript length from the other samples.
    lengthMat_notNan = replaceMissingLength(lengthMat, aveLengthSamp_g)
    lengthMat_notNan_np = lengthMat_notNan.values

    if countsFromAbundance != "no":
      countsMat_g = makeCountsFromAbundance(countsMat_g_np, abundanceMat_g_np, lengthMat_notNan_np, countsFromAbundance)
      countsMat_g = pd.DataFrame(countsMat_g, index = abundanceMat_g.index, columns = folder_names)

    if "infReps" in txi: #Comprobar con la función de Dafne (no la he encontrado).
      if varReduce:
        vars = {}
        for key, values in infReps_dict.items():
          np_i = infReps_dict[key].values
          vars_i = np.var(np_i, axis=1, ddof=1)
          vars_df = pd.DataFrame(vars_i, index = out['abundance'].index)
          vars[key] = vars_df

        vars_g = pd.concat((value for key, value in vars.items()), axis=1, keys = folder_names)


        #InfReps_vars ={} #Inicializar un diccionario para las varianzas
        #for key, values in infReps_dict.items():
          #vars_i = rowVars(infReps_dict[key])
          #infReps_vars[key] = vars_i

        out = {"abundance":abundanceMat_g,
               "counts":countsMat_g,
               "variance":vars_g,
               "length":lengthMat_notNan,
               "countsFromAbundance":countsFromAbundance}
      else:
        out = {"abundance":abundanceMat_g,
               "counts":countsMat_g,
               "infReps":infReps_gDict,
               "length":lengthMat_notNan,
               "countsFromAbundance":countsFromAbundance}
    else:
      out = {"abundance":abundanceMat_g,
             "counts":countsMat_g,
             "length":lengthMat_notNan,
             "countsFromAbundance":countsFromAbundance}

  return out