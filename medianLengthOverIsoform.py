def medianLengthOverIsoform(length, tx2gene, ignoreTxVersion = True, ignoreAfterBar = False):
  import re
  import pandas as pd

  txId = raw[txIdCol]
  length_df = pd.DataFrame(length, index=txId)
  lenColNum = len(length[0])

  if ignoreTxVersion:
    #Eliminamos todas las versiones detrás del punto para sacar las isoformas de los tránscritos.
    txId = txId.apply(lambda x: re.sub("\\..*", "", str(x)))

  elif ignoreAfterBar:
    #Eliminamos las versiones detrás de la barra para sacar las isoformas.
    txId= txId.apply(lambda x: re.sub("\\|.*", "", str(x)))

  #Crear un df de las Ids
  txId_df=pd.DataFrame({"trásncritosId":txId})

  # Cambiar el nombre de las columnas del DF de las anotaciones
  tx2gene = tx2gene.set_axis(["tx","gene"], axis='columns')

  #Creamos una serie de pandas sin isoformas de los tránscritos
  tx2gene_tx = tx2gene["tx"].apply(lambda x: re.sub("\\..*", "", str(x)))

  if not txId.isin(tx2gene_tx).all():
    #tx2gene_tx.isin(txId).all():
    print("Not all values of raw[txIdCol] column are present in 'tx' column of tx2gene.",
          "There will be transcrips with NaN gene anotation.")

  # Crear una copia del DataFrame original y asignar los valores sin isoformas de la columna 'tx'
  tx2gene_modified = tx2gene.copy()
  tx2gene_modified['tx'] = tx2gene_tx

  #Limpiamos el tx2gene de tráscritos duplicados y cambiamos sus nombres de las columas:
  tx2gene_cleaned = cleanTx2gene(tx2gene_modified)

  #Calculamos la intersección entre txId y tx2gene para sacar los tráscritos anotados presentes
  #interseccion = pd.Series(list(set(tx2gene_cleaned['tx']) and set(txId)))
  interseccion = set(tx2gene_cleaned['tx']).intersection(set(txId))
  interseccion = pd.Series(list(interseccion))

  #Creamos un df con la intersección:
  interseccion_df = pd.DataFrame({'elemento': interseccion})

  #Creamos el df con las filas con intersección: conservamos solo los genes que corresponden con interseccion_df
  merged_df = pd.merge(interseccion_df, tx2gene_cleaned, left_on= 'elemento', right_on='tx', how='left')

  # Extraer las filas con valores NaN en alguna columna
  filas_con_nan = merged_df[merged_df.isna().any(axis=1)]

  #Se limpia el df de tránscritos con anotaciones de gen == NaN (puede ser redundante)
  df_cleaned = merged_df.dropna()

  ave_length = length_df.mean(axis=1)
  ave_length_df = pd.DataFrame({'Id_transcript': ave_length.index, 'Valor': ave_length.values})

  # Agrupar los índices por el valor de la columna 'gene'
  grupos = df_cleaned.groupby('gene').groups
  longitud=len(grupos.keys())


  # Iterar sobre cada lista en el diccionario y convertir cada elemento a entero
  for clave, lista in grupos.items():
      grupos[clave] = [int(elemento) for elemento in lista]

  # Como sobran pares clave-valor, crear una lista de claves a eliminar
  claves_a_eliminar = []

  # Iterar sobre el diccionario
  for clave, valor in grupos.items():
      if not valor:
          claves_a_eliminar.append(clave)

  # Eliminar las claves del diccionario
  for clave in claves_a_eliminar:
      del grupos[clave]

  # Comprobar que hemos borrado elementos del diccionario
  longitud2=len(grupos.keys())

  indices_agrupados = list(grupos.values())
  claves_grupo = list(grupos.keys())

  print("La nueva longitud tras eliminar los genes sin índice:",longitud2)

  # Inicializar una lista para almacenar las medianas por grupo
  medianas_por_grupo = []

  # Calcular la mediana para cada grupo
  for indices_grupo in indices_agrupados:
      # Seleccionar los elementos de la Serie correspondientes al grupo actual
      elementos_grupo = ave_length_df["Valor"].iloc[indices_grupo]
      # Calcular la mediana para el grupo actual
      mediana_grupo = elementos_grupo.median()
      # Agregar la mediana a la lista de medianas
      medianas_por_grupo.append(mediana_grupo)

  # Crear un DataFrame con los datos de las medianas asociadas a su gen
  df_medianas = pd.DataFrame({'Gen': claves_grupo, 'Mediana': medianas_por_grupo})
  # Convertir la columna de medianas a un array de numpy y calculamos su longitud
  medianas_np = medianas['Mediana'].values.reshape(-1, 1)
  # Calcular la longitud del vector de las medianas
  lenMedianas_np = len(medianas_np)
  # Concatenar el array tantas veces como muestras tengamos
  medReps = np.tile(medianas_np, lenColNum)
  # Reordenar el array concatenado en una matriz de (nº genes x nº muestras)
  matMedianas = np.reshape(medReps,(lenMedianas_np, lenColNum))

  return matMedianas