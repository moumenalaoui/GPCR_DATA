**GPCR Structure-Based Feature Extraction & Machine Learning Preparation**

**Researcher:** Moumen Alaoui  
**Supervisor:** Dr. Alex Dickson, Samik Bose 
**Target Receptor:** 5-HT7 (Serotonin Receptor Subtype)  
**Goal:** Prepare molecular interaction and structural descriptors for downstream machine learning models predicting GPCR ligand behavior.

---

### Objectives
- Clean and preprocess ligand–receptor complex files.
- Generate three feature groups:
  - **INT_Feat**: Protein–ligand interaction fingerprints.
  - **LIG_Feat**: Ligand molecular descriptors.
  - **POCK_Feat**: Binding pocket characteristics.
- Merge the features into a unified matrix for ML model input.
- Attempt classification model using the open-source **GPCR-IPL_Score** framework.

---

### Tools Used
- **RDKit**: Ligand descriptor extraction.
- **PLIP**: Interaction fingerprinting.
- **Python (with Pandas, XGBoost, Tensorflow)**: Data manipulation and modeling.
- **Anaconda environments**: Environment reproducibility.

---

### Pipeline Summary
1. **Ligand Preparation**
   - Started with `5HT.smi` SMILES file.
   - Generated ligand descriptors using RDKit (e.g., molecular weight, LogP).

2. **Protein–Ligand Complex Preparation**
   - Cleaned and formatted PDB files (`XP_500_bestposes_complex.pdb`) using `fix_protein_atoms.py` and `fix_pdb_format.py`.

3. **Feature Generation**
   - INT_Feat: Extracted with `extract_interactions.py` via PLIP.
   - POCK_Feat: Pocket atom count from `extract_pocket.py`.
   - LIG_Feat: RDKit-based descriptors via `convert_5HT.py`.

4. **Feature Integration**
   - Merged all features into `features_matrix.csv`.

5. **Model Testing**
   - Attempted to run `xgb_classifier()` from the GPCR-IPL_Score repo.
   - Limited to one sample, so full training/testing was not possible.

---

### Output
- **Final Output File:** `features_matrix.csv`
- **GitHub Repository:** https://github.com/moumenalaoui/GPCR-IPL_Score-modified
- **Directory Structure:** All relevant scripts, input files, and outputs are stored in the `GPCR_DATA/` folder.

---

### Challenges
- Only one ligand–receptor sample available, preventing meaningful ML training.
- Compatibility and environment setup for Tensorflow and XGBoost required debugging.
- No pre-trained model available in GPCR-IPL_Score repo.

---

### Future Work
- Dock more ligands to 5-HT7 to generate a multi-sample dataset.
- Explore ligand clustering for dataset diversity.
- Retrain the XGBoost model using a larger `features_matrix.csv`.
- Integrate explainable ML metrics (e.g., SHAP from repo).

---

Prepared by: **Moumen Alaoui**  
Date: April 3, 2025

