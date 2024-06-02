def load_tx2gene(tx2gene_file):
  import pandas as pd
  tx2gene = pd.read_csv(tx2gene_file)
  return tx2gene

def cleanTx2gene(tx2gene):
  tx2gene = tx2gene.set_axis(["tx","gene"], axis='columns')

  if tx2gene.duplicated(subset=["tx"]).any():
    # Si hay filas duplicadas en la columna 'tx', elimínalas
    duplicados = tx2gene[tx2gene.duplicated(subset=['tx'], keep='first')]
    print("Removiendo las siguientes filas duplicadas de la columna 'tx':",
          duplicados)
    tx2gene = tx2gene.drop_duplicates(subset=['tx'], keep='first')

  tx2gene['gene'] = pd.Categorical(tx2gene['gene'])
  tx2gene['tx'] = pd.Categorical(tx2gene['tx'])

  return tx2gene