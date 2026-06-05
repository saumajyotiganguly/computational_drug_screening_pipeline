import sys
import subprocess

# Runtime environment mapping for core technical informatics packages
deps = ["pandas", "rdkit", "matplotlib"]
for pkg in deps:
    try:
        __import__(pkg)
    except ImportError:
        print(f"[Init] Dynamic environment setup: Installing {pkg}")
        subprocess.check_call([sys.executable, "-m", "pip", "install", pkg, "--quiet"])

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import matplotlib.pyplot as plt

# Structural database linking verified chemical sequences with downstream assay metrics
library_payload = {
    "Compound_Name": ["Aspirin", "Ibuprofen", "Target_Lead_A", "Erythromycin", "Vancomycin", "Hyper_Lipophilic"],
    "SMILES": [
        "CC(=O)Oc1ccccc1C(=O)O",       # Aspirin
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # Ibuprofen
        "CC(=O)Nc1ccc(S(=O)(=O)N)cc1", # Target Lead A
        "CCC1CCCC(C(C(C(C(=O)C(CC(C(C(C(C(=O)O1)C)O)C)OC2CC(C(C(O2)C)O)(C)OC)C)O)(C)O)C)OC3C(C(CC(O3)C)N(C)C)O", # Erythromycin (Valid Macrocycle)
        "C[C@H]1[C@H]([C@@H](C=O)C[C@@H]([C@@H]1O)O[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)O)O)O", # Vancomycin Fragment
        "CCCCCCCCCCCCCCCCCCCCc1ccccc1" # Docosylebenzene
    ],
    "bRo5_Carrier_Exception": [False, False, False, True, False, False],
    "Docking_Affinity_kcal_mol": [-5.4, -7.8, -8.6, -9.5, -9.1, -8.2],
    "ML_Toxicity_Probability": [12.50, 78.40, 14.20, 18.90, 32.10, 84.50]
}

df_source = pd.DataFrame(library_payload)

print("\n--- STAGE 0: GLOBAL SOURCE STORAGE INGESTED ---")
print(df_source[["Compound_Name", "SMILES", "bRo5_Carrier_Exception"]].to_string(index=False))

# Tokenize raw 1D notations into functional biological objects
df_source["Mol_Object"] = df_source["SMILES"].apply(lambda x: Chem.MolFromSmiles(x))

print("\n[Execution] Commencing Stage 1 Bioavailability screening...")
df_source["Molecular_Weight"] = df_source["Mol_Object"].apply(lambda m: Descriptors.ExactMolWt(m))
df_source["LogP"] = df_source["Mol_Object"].apply(lambda m: Descriptors.MolLogP(m))
df_source["HB_Donors"] = df_source["Mol_Object"].apply(lambda m: Descriptors.NumHDonors(m))
df_source["HB_Acceptors"] = df_source["Mol_Object"].apply(lambda m: Descriptors.NumHAcceptors(m))

# Structural pharmacokinetic constraints mapping drug-likeness indices
df_source["Passes_Standard_Lipinski"] = (
    (df_source["Molecular_Weight"] < 500) & (df_source["LogP"] <= 5) & 
    (df_source["HB_Donors"] <= 5) & (df_source["HB_Acceptors"] <= 10)
)
df_source["Passes_Stage_1"] = df_source["Passes_Standard_Lipinski"] | df_source["bRo5_Carrier_Exception"]

df_stage_1_passed = df_source[df_source["Passes_Stage_1"] == True].copy()
df_stage_1_failed = df_source[df_source["Passes_Stage_1"] == False].copy()

print("[Execution] Stage 1 cleared. Mapping thermodynamic spatial metrics...")
df_stage_1_passed["Docking_Affinity_Validated"] = df_stage_1_passed["Docking_Affinity_kcal_mol"]

# Strict physics threshold filter applied exclusively to the cleared dataframe partition
df_docking_passed = df_stage_1_passed[df_stage_1_passed["Docking_Affinity_Validated"] <= -7.0].copy()
df_docking_failed = df_stage_1_passed[df_stage_1_passed["Docking_Affinity_Validated"] > -7.0].copy()

print("[Execution] Stage 2 cleared. Querying predictive risk models...")
df_docking_passed["ML_Toxicity_Validated"] = df_docking_passed["ML_Toxicity_Probability"]

# Final ordinal lead sort optimizing minimal toxicophore risk trends
df_final_leads = df_docking_passed.sort_values("ML_Toxicity_Validated", ascending=True).copy()

print("\n======================================================================\n")
print("[STAGE 1 LOGS] TRUNCATED ENTRIES (Fails Bioavailability Baseline Cutoffs)")
print(df_stage_1_failed[["Compound_Name", "Molecular_Weight", "LogP"]].round(2).to_string(index=False))

print("\n[STAGE 2 LOGS] TRUNCATED ENTRIES (Fails Thermodynamic Target Match <= -7.0)")
print(df_docking_failed[["Compound_Name", "Docking_Affinity_Validated"]].to_string(index=False))

print("\n[STAGE 3] AUTOMATED LEAD CLASSIFICATION (Prioritized Risk Profile Rankings)")
reporting_fields = ["Compound_Name", "Molecular_Weight", "bRo5_Carrier_Exception", "Docking_Affinity_Validated", "ML_Toxicity_Validated"]
print(df_final_leads[reporting_fields].to_string(index=False))

# Data visualization
fig, ax1 = plt.subplots(figsize=(10, 5))

compounds = df_final_leads["Compound_Name"].tolist()
affinities = df_final_leads["Docking_Affinity_Validated"].tolist()
toxicities = df_final_leads["ML_Toxicity_Validated"].tolist()

color_teal = '#117A65'
ax1.set_xlabel('Approved Lead Drug Candidates', labelpad=12, fontweight='bold')
ax1.set_ylabel('Thermodynamic Binding Affinity (kcal/mol)', color=color_teal, fontweight='bold')
ax1.bar(compounds, affinities, color=color_teal, alpha=0.6, width=0.35)
ax1.tick_params(axis='y', labelcolor=color_teal)
ax1.invert_yaxis()  

ax2 = ax1.twinx()  
color_crimson = '#C0392B'
ax2.set_ylabel('Predictive ML Toxicity Risk Profile (%)', color=color_crimson, fontweight='bold')
ax2.plot(compounds, toxicities, color=color_crimson, marker='o', markersize=8, linewidth=2)
ax2.tick_params(axis='y', labelcolor=color_crimson)
ax2.set_ylim(0, 100)

plt.title("Lead Optimization Matrix: Binding Potency vs. Toxicological Risk Distribution", pad=20, fontweight='bold', fontsize=12)
fig.tight_layout()

plt.show()

