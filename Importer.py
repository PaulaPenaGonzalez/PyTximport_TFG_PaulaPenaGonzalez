def importer(x): #Definimos una función para importar los datos, donde X es cada archivo en files
  #Los tipos de objetos que esperamos encontrar en cada columna
  column_types = {"col1": str, "col2": int, "col3": float, "col4": float, "col5": float}
  data = pd.read_csv(x, sep="\t", dtype=column_types) #Leer los archivos con los tipos de columna especificados
  return data


def readInfRepFish(
        fish_dir):  # Definimos la función readInfRepFish para importar las réplicas inferenciales de SailFish.
  import os
  import numpy as np
  import gzip
  import sys
  import json
  aux_dir = "aux_info"
  cmd_info_jsonPath = os.path.join(fish_dir, "cmd_info.json")  # Fish_dir es el directorio de salida de SailFish

  if not os.path.exists(cmd_info_jsonPath):  # Comprobar que el archivo 'cmd_info.json' existe.
    # El archivo no existe imprime un error
    print("The file cmd_info.json doesn't exist in the specified path.")
  else:
    # El archivo existe, puedes continuar con el código
    print("The file cmd_info.json exists in the specified path.")

    # Abrir el archivo como diccionario y comprobar que la key 'AuxDir' existe. ¿Dentro de else?
    with open(cmd_info_jsonPath, 'r') as cmd_info:
      dict_cmd_info_json = json.load(cmd_info)  # Abrimos cmd_info como un diccionario
      aux_dir = dict_cmd_info_json[
        "auxDir"]  # Accedemos a la key del diccionario y guardamos su valor en una variable llamada 'aux_Dir

  # Concatenar el directorio de entrada con el valor de la key del diccionario (que es "aux_info")
  aux_path = os.path.join(fish_dir, aux_dir)
  if not os.path.exists(aux_path):
    print("The folder 'aux_info' doesn't exist in the specified path.")
  else:
    print("The folder 'aux_info' exists in the specified path.")
    # Crear ruta al archivo 'meta_info.json' dentro de la carpeta 'aux_info'
    meta_info_jsonPath = os.path.join(aux_path, "meta_info.json")
    # Carga los datos JSON en un diccionario de Python
    with open(meta_info_jsonPath, 'r') as meta_info:
      dict_meta_info = json.load(meta_info)
      # Comprobar las versiones de Salmon y Sailfish
      if "salmon_version" in dict_meta_info:
        salmon_version = dict_meta_info["salmon_version"]
        print("Searching for compatible salmon version")
        if salmon_version >= "0.8.0":
          print(f"Salmon version:{salmon_version}, is compatible")
        else:
          sys.exit("Salmon version >=0.8.0 is required")
      else:
        print("Salmon version is not accesible")

      if "sailfish_version" in dict_meta_info:
        sailfish_version = dict_meta_info["sailfish_version"]
        print("Searching for compatible sailfish version")
        if sailfish_version >= "0.9.0":
          print(f"Sailfish version:{sailfish_version}, is compatible")
        else:
          sys.exit("Sailfish version >=0.9.0 is required")
      else:
        print("Sailfish version is not accesible")

    sampType = None
    # Verificar si se ha registrado explícitamente el tipo de muestra posterior.
    if "samp_type" in dict_meta_info:
      sampType = dict_meta_info["samp_type"]
      print(f"Sample type is:{sampType}")

    # Cargar datos de arranque
    knownSampleTypes = ["gibbs", "bootstrap"]
    numBoot = dict_meta_info["num_bootstraps"]

    if numBoot > 0:
      if "num_valid_targets" in dict_meta_info:
        dict_meta_info["num_targets"] = dict_meta_info["num_valid_targets"]

      # Calcular cantidad total de datos de arranque que se esperan tener en el conjunto de datos
      expected_n = dict_meta_info["num_targets"] * dict_meta_info["num_bootstraps"]

      try:
        try:
          # Estaría bien crear una función con la apertura del archivo binario para no repetir código
          with gzip.open(os.path.join(aux_path, "bootstrap", 'bootstraps.gz'), "rb") as boot_Con:
            dt = np.dtype('float64')  # Tipo de dato al que queremos convertir el archivo binario
            dt = dt.newbyteorder(
              '<')  # little-endian (bit menos significativo primero):< y big-endian (más significativo primero):>
            boots_In = np.frombuffer(boot_Con.read(), dtype=dt, count=expected_n)
            tamaño = boots_In.size
            if tamaño != expected_n:
              raise ValueError("The number of elements read was not the expected one")
            else:
              print(
                "'Bootstraps.gz' was unzipped and saved successfully in 'boots' variable as an array of double precision floats")
              boots = boots_In.reshape(dict_meta_info["num_targets"], dict_meta_info["num_bootstraps"], order='F')
              vars = np.var(boots, axis=1, ddof=1)
              vars_boots_dict = {"vars": vars, "reps": boots}
              return vars_boots_dict

        except:
          with gzip.open(os.path.join(aux_path, "bootstrap", 'bootstraps.gz'), "rb") as boot_Con:
            boots_In = np.frombuffer(boot_Con.read(), dtype=int, count=expected_n)
            tamaño = boots_In.size
            if tamaño != expected_n:
              raise ValueError("The number of elements read was not the expected one")
            else:
              # Crear una matriz de tamaño de nxn: tráncritos x bootstraps
              boots = boots_In.reshape(dict_meta_info["num_targets"], dict_meta_info["num_bootstraps"], order='F')
              vars = np.var(boots, axis=1, ddof=1)  # Calcula la varianza por filas de la matriz
              vars_boots_dict = {"vars": vars, "reps": boots}  # Guardar varianzas y matriz en un diccionario
              return vars_boots_dict
              print("'Bootstraps.gz' was unzipped and saved successfully in 'boots' variable as an array of integers")


      except ValueError as ve:
        print("ValueError:", ve)

      except FileNotFoundError:
        # Manejar el caso en que el archivo no existe
        boots_In = None
        print("File 'bootstraps.gz' not found")

      except TypeError as te:
        print("TypeError:", te)

      except Exception as e:
        # Cualquier otra excepción
        print("Error: The file could not be translated correctly to float or integer", e)

    else:
      print("There are no Bootstraps: variances could not be calculated")
      vars_boots_dict = None
      return vars_boots_dict


def read_kallito_h5(fpath):
  import h5py
  import pandas as pd
  import numpy as np
  with h5py.File(fpath,'r') as archivo:
    todas_las_claves = []
    archivo.visit(todas_las_claves.append)

    #Si existe un grupo llamado 'est_counts' llenar un array de 1D con sus valores
    if 'est_counts' in todas_las_claves:
      numCounts=len(archivo['est_counts'])
      counts= np.zeros((numCounts))
      counts[:] = archivo['est_counts']
    else:
      raise ValueError("The group 'est_counts' was not found in abundance.h5")

    #Si existe un grupo llamado 'aux/ids' llenar un array de 1D con sus valores
    if 'aux/ids' in todas_las_claves:
      numIds=len(archivo['aux']['ids'])
      ids= pd.Series(archivo['aux']['ids'])
      ids=ids.apply(lambda x: x.decode('utf-8'))
    else:
      raise ValueError("The group 'aux/ids' was not found in abundance.h5")

    if 'aux/eff_lengths' in todas_las_claves:
      numEfflens=len(archivo['aux']['eff_lengths'])
      efflens = np.zeros((numEfflens))
      efflens[:]=archivo['aux']['eff_lengths']
    else:
      raise ValueError("The group 'aux/eff_lengths' was not found in abundance.h5")

    if not len(efflens)== len(ids) and len(counts) == len(ids):
      raise ValueError()

    result = pd.DataFrame({
    'target_id': ids,
    'eff_length': efflens,
    'est_counts': counts,
    'stringsAsFactors' : False})

    normfac = 1e6/(sum(counts/efflens))

    result['tpm'] = normfac * (result['est_counts'] / result['eff_length'])
    return result


def readInfRepKallisto(bearDir):
  import h5py
  h5file = os.path.join(bearDir, 'abundance.h5')
  if os.path.exists(h5file):
    print(f"The file {h5file} exists")
  else:
    return None
  with h5py.File(h5file,'r') as archivo:
    grupos=list(archivo.keys())
    if 'bootstrap'in grupos:
      print("El grupo 'bootstraps' si que existe en la lista grupos")
      grupos_bootstraps = archivo['bootstrap']
      numBoot = len(grupos_bootstraps)
      print(numBoot)
      if numBoot > 0:
        #contenido_grupo_bootstrap = grupos_bootstraps[:]
        #print(contenido_grupo_bootstrap)
        bs0=archivo['bootstrap']['bs0']
        print(bs0)
        numTx=len(archivo['bootstrap']['bs0'])
        print(numTx)
        bootMat= np.zeros((numTx, numBoot))
        print(bootMat.shape)
        for n in range(numBoot):
          column = archivo['bootstrap'][f'bs{n}']
          bootMat[:,n]= column
        print(bootMat)
        vars = np.var(bootMat, axis=1)
        print(vars)
        vars_boots_dict = {"vars": vars, "reps": bootMat}
        print(vars_boots_dict)
        return vars_boot_dict
      else:
        return None

