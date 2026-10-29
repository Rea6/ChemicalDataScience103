#- Prepare the data from ChEMBL: Load the csv as DataFrame and select only the following columns: 'Molecule ChEMBL ID', 'Smiles', 'Standard Type', 'Standard Relation', 'Standard Value', 'Standard Units', 'Assay ChEMBL ID'.
#- Rename the columns (e.g. avoid spaces, more convenient names). 
#- You can drop all rows which contain no Smiles code (make sure to comment!).
#- Select from the DataFrame entries with a specific molecular weight range (e.g. 200 to 400). You can change this value later to see how this might affect your hypothesis.
#- Select a couple of hundred entries (pseudo)randomly from the entire DataFrame (don't use too many, e.g. >1000, otherwise, the pubchem query may take very long). You can use a pandas method for that: 
 #  ```python
  # df_sample = df.sample(n=300, random_state=42)
   #```
#- Export this DataFrame to a .csv file.



import pandas as pd
from tqdm import tqdm
import time
from pathlib import Path

df = pd.read_csv(Path(__file__).parent / "CHEMBL_task12.csv", sep=";")

df2= pd.DataFrame()
df2['MolID']= df['Molecule ChEMBL ID']
df2['Smiles']=df['Smiles']
df2['Std_Type']=df['Standard Type']
df2['Std_Relation']=df['Standard Relation']
df2['Std_Value']=df['Standard Value']
df2['Std_Units']=df['Standard Units']
df2['Assay_ID']=df['Assay ChEMBL ID']

df2.describe()
