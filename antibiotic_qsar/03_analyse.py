"""Merge, standardize, visualize and test MIC–XLogP relationships.
Results are exploratory; repeated assay observations are not automatically independent.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from rdkit import Chem
from rdkit.Chem import Draw

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data'
OUT = ROOT / 'results'
OUT.mkdir(exist_ok=True)
chembl = pd.read_csv(DATA / 'chembl_sample.csv')
pubchem = pd.read_csv(DATA / 'antibio_pubchem_descr.csv')
print('ChEMBL:', chembl.shape, 'PubChem:', pubchem.shape)
print('PubChem duplicated SMILES:', pubchem['smiles'].duplicated().sum())
# Avoid accidental many-to-many expansion caused by repeated PubChem SMILES.
pubchem = pubchem.drop_duplicates(subset=['smiles'], keep='first')
df = chembl.merge(pubchem, on='smiles', how='inner', validate='many_to_one')
print('Merged observations:', len(df))
print('Original MIC units:\n', df['mic_units'].value_counts(dropna=False))

# Convert exact concentration measurements to a common unit: micrograms per litre.
# mg/L = ug/mL = 1000 ug/L; mg/mL = 1,000,000 ug/L.
# molar concentrations require molecular weight (g/mol): 1 uM = MW ug/L.
# Unknown units are deliberately discarded, not treated as equivalent.
def convert_mic(row):
    unit = str(row['mic_units']).strip().lower().replace('µ', 'u').replace('μ', 'u').replace(' ', '')
    value = row['mic_value']
    weight = row['mol_weight_rdkit']
    scale = {
        'ug/l': 1., 'mcg/l': 1., 'ng/ml': 1.,
        'mg/l': 1000., 'ug/ml': 1000., 'mcg/ml': 1000.,
        'mg/ml': 1_000_000., 'ng/l': 0.001,
        'nm': weight * 0.001, 'um': weight, 'mm': weight * 1000.,
        'm': weight * 1_000_000.,
    }
    return float(value) * scale[unit] if unit in scale else np.nan

df['mic_ug_l'] = df.apply(convert_mic, axis=1)
print('Unconvertible rows:', df['mic_ug_l'].isna().sum())
# Exclude missing / nonpositive MIC and descriptors, and ambiguous units.
for col in ['XLogP', 'TPSA', 'MolecularWeight', 'mic_ug_l']:
    df[col] = pd.to_numeric(df[col], errors='coerce')
df = df.dropna(subset=['smiles', 'XLogP', 'mic_ug_l']).copy()
df = df[(df['mic_ug_l'] > 0) & np.isfinite(df['mic_ug_l'])].copy()
# Use one median assay measurement per compound, avoiding pseudoreplication in the t-test.
# This simplified aggregation loses assay-specific variation; retain raw merged data for audit.
df.to_csv(OUT / 'merged_clean_measurements.csv', index=False)
df = (df.groupby('smiles', as_index=False)
      .agg(chembl_id=('chembl_id', 'first'), mic_ug_l=('mic_ug_l', 'median'),
           XLogP=('XLogP', 'first'), TPSA=('TPSA', 'first'),
           MolecularWeight=('MolecularWeight', 'first'),
           assay_count=('assay_id', 'nunique')))
df['log10_mic'] = np.log10(df['mic_ug_l'])
df.to_csv(OUT / 'cleaned_compounds.csv', index=False)
print('Independent compounds:', len(df))
print(df.describe(include='all').to_string())
print('Null values:\n', df.isna().sum())
print('10 smallest MIC values:\n', df.nsmallest(10, 'mic_ug_l')[['chembl_id', 'mic_ug_l', 'XLogP']])

# Top ten compounds: mol structures labelled with MIC concentration.
best = df.nsmallest(10, 'mic_ug_l')
molecules, legends = [], []
for _, row in best.iterrows():
    mol = Chem.MolFromSmiles(row['smiles'])
    if mol is not None:
        molecules.append(mol)
        legends.append(f"{row['chembl_id']} | {row['mic_ug_l']:.3g} ug/L")
if molecules:
    Draw.MolsToGridImage(molecules, molsPerRow=5, subImgSize=(270, 230),
                         legends=legends).save(str(OUT / 'top_10_molecules.png'))

# Simple EDA distributions and scatterplot; MIC is log transformed to handle large ranges.
for column, filename in [('XLogP', 'xlogp_histogram.png'), ('log10_mic', 'log_mic_histogram.png')]:
    plt.figure(figsize=(7, 5))
    plt.hist(df[column].dropna(), bins=25)
    plt.xlabel(column)
    plt.ylabel('Number of molecules')
    plt.tight_layout()
    plt.savefig(OUT / filename, dpi=160)
    plt.close()

if len(df) < 4:
    print('Too few valid compounds for reliable t-tests.')
else:
    # Exploratory split at median XLogP; midpoint avoids arbitrarily setting a cutoff.
    # A cutoff tuned to maximize significance would invalidate a naive p-value.
    cutoff = df['XLogP'].median()
    print('Exploratory XLogP threshold (median):', cutoff)
    plt.figure(figsize=(8, 5))
    plt.scatter(df['XLogP'], df['log10_mic'], alpha=0.7)
    plt.axvline(cutoff, color='black', linestyle='--', label=f'Median XLogP = {cutoff:.2f}')
    plt.xlabel('PubChem XLogP')
    plt.ylabel('log10(MIC in ug/L)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / 'xlogp_vs_log_mic.png', dpi=160)
    plt.close()

    def test_groups(frame, label):
        lower = frame.loc[frame['XLogP'] <= cutoff, 'log10_mic'].dropna()
        upper = frame.loc[frame['XLogP'] > cutoff, 'log10_mic'].dropna()
        print(f'\n{label}: Low-XLogP n={len(lower)}, high-XLogP n={len(upper)}')
        if len(lower) < 2 or len(upper) < 2:
            print('Not enough observations in both groups to perform Welch t-test.')
            return
        result = ttest_ind(lower, upper, equal_var=False)
        print(f'Mean log10(MIC): low={lower.mean():.3f}, high={upper.mean():.3f}')
        print(f'Welch t={result.statistic:.4f}, p={result.pvalue:.5g}')
        if result.pvalue < .05:
            print('At alpha=0.05, the observed mean difference is statistically significant.')
        else:
            print('At alpha=0.05, insufficient evidence for a mean difference.')
        print('A lower group mean log10(MIC) corresponds to a lower MIC.')

    test_groups(df, 'All valid compounds')
    test_groups(df[df['mic_ug_l'] < 10], 'MIC < 10 ug/L (as assigned)')
    print('Important: Selecting on the MIC outcome changes the population and may bias inference.')
