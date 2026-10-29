Session 12 – E. coli MIC / molecular descriptors

1. Export ChEMBL target CHEMBL354, MIC activity, with molecular-weight filter up to 600.
2. Save the downloaded export as: data/chembl_ecoli_mic.csv
3. Install packages in your Python environment:
   pip install pandas numpy scipy matplotlib rdkit pubchempy tqdm
4. Run scripts from this folder in this order:
   python 01_prepare_chembl.py
   python 02_pubchem_search.py
   python 03_analyse.py

Outputs:
 data/chembl_sample.csv
 data/antibio_pubchem_descr.csv
 results/cleaned_compounds.csv
 results/merged_clean_measurements.csv
 results/top_10_molecules.png
 results/xlogp_vs_log_mic.png
 results/xlogp_histogram.png
 results/log_mic_histogram.png

Notes:
 - The ChEMBL CSV has not been provided, so the workflow is not yet executed on real data.
 - 200 to 400 g/mol sampling uses RDKit molecular weights, prior to PubChem.
 - Non-equality or missing MIC measurements, unresolved compounds, and unknown units are omitted.
 - Measurements of the same SMILES are summarized using the median.
 - The second test uses literal MIC < 10 ug/L per the assignment; confirm that
   the instructor did not mean 10 ug/mL.
 - Assay conditions/strain differences are not controlled, so results are exploratory.
 - If the experiment uses one threshold chosen from the plotted data, p-values
   must be interpreted as exploratory rather than confirmatory.
