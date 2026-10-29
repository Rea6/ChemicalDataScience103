"""Adaptation of the provided pubchem_search.py; requires internet access."""
from pathlib import Path
import time
import pandas as pd
import pubchempy as pcp
from tqdm import tqdm

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data'
chembl_df = pd.read_csv(DATA / 'chembl_sample.csv')
# Query unique SMILES only; the same SMILES can appear for several ChEMBL assays.
smiles_list = chembl_df['smiles'].dropna().astype(str).drop_duplicates().tolist()
records = []
for smi in tqdm(smiles_list, desc='PubChem lookup'):
    try:
        compounds = pcp.get_compounds(smi, namespace='smiles')
        if not compounds:
            continue
        c = compounds[0]
        records.append({
            'smiles': smi,
            'CID': c.cid,
            'MolecularWeight': c.molecular_weight,
            'XLogP': c.xlogp,
            'TPSA': c.tpsa,
            'HBondDonorCount': c.h_bond_donor_count,
            'HBondAcceptorCount': c.h_bond_acceptor_count,
            'RotatableBondCount': c.rotatable_bond_count,
            'HeavyAtomCount': c.heavy_atom_count,
        })
        time.sleep(0.2)  # Pause briefly between PubChem API requests.
    except Exception as exc:
        print(f'Failed to find {smi[:35]}...: {exc}')

out = DATA / 'antibio_pubchem_descr.csv'
pd.DataFrame(records, columns=['smiles','CID','MolecularWeight','XLogP','TPSA',
    'HBondDonorCount','HBondAcceptorCount','RotatableBondCount','HeavyAtomCount']).to_csv(out, index=False)
print(f'Saved {len(records)} descriptor records to {out}')
