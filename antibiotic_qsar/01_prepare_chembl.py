"""Prepare a reproducible 300-row ChEMBL sample for PubChem lookups.
Place the exported ChEMBL CSV in data/chembl_ecoli_mic.csv.
"""
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data'
DATA.mkdir(exist_ok=True)
INPUT = DATA / 'chembl_ecoli_mic.csv'
OUTPUT = DATA / 'chembl_sample.csv'

# sep=None automatically detects comma- or semicolon-delimited CSV exports.
df = pd.read_csv(INPUT, sep=None, engine='python')
required = ['Molecule ChEMBL ID', 'Smiles', 'Standard Type', 'Standard Relation', 'Standard Value', 'Standard Units', 'Assay ChEMBL ID']
missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f'Columns not found: {missing}. Available: {list(df.columns)}')

# Keep only the seven columns required in the assignment.
df = df[required].rename(columns={
    'Molecule ChEMBL ID': 'chembl_id', 'Smiles': 'smiles',
    'Standard Type': 'standard_type', 'Standard Relation': 'relation',
    'Standard Value': 'mic_value', 'Standard Units': 'mic_units',
    'Assay ChEMBL ID': 'assay_id'
})
print('Imported rows:', len(df))
# Remove entries without SMILES, since descriptors cannot be calculated for them.
df = df.dropna(subset=['smiles']).copy()
df = df[df['smiles'].astype(str).str.strip().ne('')].copy()
df['mic_value'] = pd.to_numeric(df['mic_value'], errors='coerce')
# Keep only quantitative MIC observations with exact ('=') rather than censored ('>'/'<') results.
df['relation'] = df['relation'].fillna('=').astype(str).str.strip()
df = df[(df['standard_type'].astype(str).str.upper() == 'MIC') &
        (df['relation'].isin(['=', ''])) & (df['mic_value'] > 0)].copy()

def molecular_weight(smiles):
    mol = Chem.MolFromSmiles(str(smiles))
    return Descriptors.MolWt(mol) if mol is not None else float('nan')

# Compute MW directly from SMILES because the requested 7-column selection excludes MW.
df['mol_weight_rdkit'] = df['smiles'].map(molecular_weight)
df = df[df['mol_weight_rdkit'].between(200, 400)].copy()
print('Rows with exact positive MIC and MW 200-400:', len(df))
if df.empty:
    raise ValueError('No rows remain; check SMILES, units and the molecular-weight range.')
# Reproducible random sample; multiple assay measurements may intentionally remain.
sample = df.sample(n=min(300, len(df)), random_state=42)
sample.to_csv(OUTPUT, index=False)
print(f'Saved {len(sample)} rows to {OUTPUT}')
