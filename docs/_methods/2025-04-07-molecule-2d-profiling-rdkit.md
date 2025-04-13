---
layout: post
title: "How to Standardize and Profile Chemical Structures with RDKit"
description: "Learn how to standardize and analyze chemical structures using RDKit's Python API. Step-by-step tutorial on molecule processing, descriptor calculation, and more for cheminformatics workflows."
author: Brad Dallin
catagories: methods
---

Chemical structure standardization and descriptor calculation are foundational steps in modern drug discovery pipelines and computational chemistry workflows. In this guide, I'll demonstrate essential techniques for processing, standardizing, and calculating molecular descriptors using RDKit—a powerful open-source cheminformatics toolkit widely adopted in pharmaceutical research and computational chemistry.

Computational chemists, medicinal chemists, and cheminformaticians frequently face challenges when handling molecular data from diverse sources. Variations in representation, salt forms, ionization states, and structural encoding can significantly impact downstream analyses such as QSAR modeling, virtual screening, and molecular docking. This tutorial addresses these challenges by providing a systematic approach to chemical data standardization—ensuring consistency and reproducibility across molecular datasets.

Building on my [previous post](https://brad-dallin.github.io/laptopchemistry/methods/2025-03-29-query-chembl-approved-small-molecule-drugs.html) about ChEMBL data retrieval, I'll show you how to efficiently process SMILES strings, implement a molecular preparation techniques, and calculate molecular descriptors critical for structure-property relationship studies and drug-likeness assessment. You'll learn practical RDKit functions for handling salt stripping, charge neutralization, and molecular feature extraction that can be integrated into Python-based cheminformatics pipelines.

This tutorial assumes basic familiarity with Python programming and chemical structure representation, making it ideal for computational chemists, pharmaceutical researchers, and data scientists working with chemical datasets.

I’ll be using a RDKit python environment I previously [created](https://brad-dallin.github.io/laptopchemistry/methods/2025-03-29-query-chembl-approved-small-molecule-drugs.html#python-setup) or see the [RDKit documentation](https://www.rdkit.org/docs/Install.html) for more details. The jupyter notebook for this methods post can be found [here](https://github.com/brad-dallin/laptopchemistry/blob/main/notebooks/rdkit_molecule_2d_profiling.ipynb).


### **1. Import modules**

Import `os`, `pandas`, and `rdkit` modules. `Descriptors` for calculating molecular descriptors and `rdMolStandardize` for molecule standardization are imported from `rdkit.Chem`.

Pandas display settings are adjusted to show all columns, and RDKit logging is suppressed to reduce excessive output. Finally, the versions of Pandas and RDKit are printed.

```python
# Import modules
import os
import pandas as pd
import rdkit
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

# Expand to see all columns
pd.set_option('display.max_columns', None)

# Suppress RDKit Output
rdkit.RDLogger.DisableLog('rdApp.*')

# Print versions
print(f"Pandas Version: {pd.__version__}")
print(f"RDKit Version: {rdkit.__version__}")
```
    Pandas Version: 2.2.3
    RDKit Version: 2024.09.6


### **2. Create ACS 1996 drawing function**

This is my personal preference for viewing 2D chemical structures and completely optional.

This function creates a publication-quality 2D visualization of a chemical structure using the American Chemical Society's 1996 style guidelines. It processes an RDKit molecule object by computing optimal 2D coordinates, applying professional styling parameters, and returning a ready-to-display image. The result is a standardized chemical structure representation.

```python
# Draw molecules in ACS1996 format
def show_acs1996(mol, legend=""):
    from PIL import Image
    from io import BytesIO
    from rdkit.Chem import Draw
    rdkit.Chem.rdDepictor.Compute2DCoords(mol)
    rdkit.Chem.rdDepictor.StraightenDepiction(mol)
    d2d = Draw.MolDraw2DCairo(-1,-1)
    Draw.DrawMoleculeACS1996(d2d, mol,legend=legend)
    bio = BytesIO(d2d.GetDrawingText())
    return Image.open(bio)
```


### **3. Load dataset**
Load the CSV file as a Pandas dataframe. The dataset used can be found in the notebooks folder of this repository [here](https://github.com/brad-dallin/laptopchemistry/blob/main/notebooks/data/chembl_approved_small_molecule_drugs.csv).

```python
# Paths
path_data = os.path.realpath("../data")
input_file = "chembl_approved_small_molecule_drugs.csv"
input_path = os.path.join(path_data, input_file)

# Load CSV
df = pd.read_csv(input_path)
print(df.shape)
df.head()
```
    (3517, 53)

<div style="overflow-x: auto; width: 100%;">
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }
    .dataframe tbody tr th {
        vertical-align: top;
    }
    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>chembl_molecule_id</th>
      <th>smiles</th>
      <th>molregno</th>
      <th>pref_name</th>
      <th>max_phase</th>
      <th>therapeutic_flag</th>
      <th>dosed_ingredient</th>
      <th>structure_type</th>
      <th>chebi_par_id</th>
      <th>molecule_type</th>
      <th>first_approval</th>
      <th>oral</th>
      <th>parenteral</th>
      <th>topical</th>
      <th>black_box_warning</th>
      <th>natural_product</th>
      <th>first_in_class</th>
      <th>chirality</th>
      <th>prodrug</th>
      <th>inorganic_flag</th>
      <th>usan_year</th>
      <th>availability_type</th>
      <th>usan_stem</th>
      <th>polymer_flag</th>
      <th>usan_substem</th>
      <th>usan_stem_definition</th>
      <th>indication_class</th>
      <th>withdrawn_flag</th>
      <th>chemical_probe</th>
      <th>orphan</th>
      <th>mw_freebase</th>
      <th>alogp</th>
      <th>hba</th>
      <th>hbd</th>
      <th>psa</th>
      <th>rtb</th>
      <th>ro3_pass</th>
      <th>num_ro5_violations</th>
      <th>cx_most_apka</th>
      <th>cx_most_bpka</th>
      <th>cx_logp</th>
      <th>cx_logd</th>
      <th>molecular_species</th>
      <th>full_mwt</th>
      <th>aromatic_rings</th>
      <th>heavy_atoms</th>
      <th>qed_weighted</th>
      <th>mw_monoisotopic</th>
      <th>full_molformula</th>
      <th>hba_lipinski</th>
      <th>hbd_lipinski</th>
      <th>num_lipinski_ro5_violations</th>
      <th>np_likeness_score</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>CHEMBL1200542</td>
      <td>CC(=O)OCC(=O)[C@H]1CC[C@H]2[C@@H]3CCC4=CC(=O)C...</td>
      <td>674493</td>
      <td>DESOXYCORTICOSTERONE ACETATE</td>
      <td>4.0</td>
      <td>1</td>
      <td>1</td>
      <td>MOL</td>
      <td>34671.0</td>
      <td>Small molecule</td>
      <td>1939.0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>NaN</td>
      <td>0.0</td>
      <td>-cort-; -ster-; -terone</td>
      <td>0</td>
      <td>-cort-; -ster-; -terone</td>
      <td>cortisone derivatives; steroids (androgens, an...</td>
      <td>Adrenocortical Steroid (salt-regulating)</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>372.51</td>
      <td>4.27</td>
      <td>4.0</td>
      <td>0.0</td>
      <td>60.44</td>
      <td>3.0</td>
      <td>N</td>
      <td>0.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>3.77</td>
      <td>3.77</td>
      <td>NEUTRAL</td>
      <td>372.51</td>
      <td>0.0</td>
      <td>27.0</td>
      <td>0.69</td>
      <td>372.2301</td>
      <td>C23H32O4</td>
      <td>4.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.96</td>
    </tr>
    <tr>
      <th>1</th>
      <td>CHEMBL1200728</td>
      <td>Cl.N=C(N)N</td>
      <td>674679</td>
      <td>GUANIDINE HYDROCHLORIDE</td>
      <td>4.0</td>
      <td>1</td>
      <td>1</td>
      <td>MOL</td>
      <td>32735.0</td>
      <td>Small molecule</td>
      <td>1939.0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>NaN</td>
      <td>1.0</td>
      <td>guan-</td>
      <td>0</td>
      <td>guan-</td>
      <td>antihypertensives (guanidine derivatives)</td>
      <td>NaN</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>59.07</td>
      <td>-1.16</td>
      <td>1.0</td>
      <td>3.0</td>
      <td>75.89</td>
      <td>0.0</td>
      <td>N</td>
      <td>0.0</td>
      <td>NaN</td>
      <td>12.55</td>
      <td>-1.24</td>
      <td>-3.65</td>
      <td>BASE</td>
      <td>95.53</td>
      <td>0.0</td>
      <td>4.0</td>
      <td>0.24</td>
      <td>59.0483</td>
      <td>CH6ClN3</td>
      <td>3.0</td>
      <td>5.0</td>
      <td>0.0</td>
      <td>0.32</td>
    </tr>
    <tr>
      <th>2</th>
      <td>CHEMBL1200982</td>
      <td>CCC(C)C1(CC)C(=O)[N-]C(=O)NC1=O.[Na+]</td>
      <td>674933</td>
      <td>BUTABARBITAL SODIUM</td>
      <td>4.0</td>
      <td>1</td>
      <td>1</td>
      <td>MOL</td>
      <td>NaN</td>
      <td>Small molecule</td>
      <td>1939.0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>NaN</td>
      <td>1.0</td>
      <td>-barb-</td>
      <td>0</td>
      <td>-barb-</td>
      <td>barbituric acid derivatives</td>
      <td>Sedative-Hypnotic</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>212.25</td>
      <td>0.79</td>
      <td>3.0</td>
      <td>2.0</td>
      <td>75.27</td>
      <td>3.0</td>
      <td>N</td>
      <td>0.0</td>
      <td>7.48</td>
      <td>NaN</td>
      <td>1.45</td>
      <td>1.19</td>
      <td>NEUTRAL</td>
      <td>234.23</td>
      <td>0.0</td>
      <td>15.0</td>
      <td>0.68</td>
      <td>212.1161</td>
      <td>C10H15N2NaO3</td>
      <td>5.0</td>
      <td>2.0</td>
      <td>0.0</td>
      <td>0.32</td>
    </tr>
    <tr>
      <th>3</th>
      <td>CHEMBL3989520</td>
      <td>NCCc1c[nH]cn1.O=P(O)(O)O.O=P(O)(O)O</td>
      <td>2197391</td>
      <td>HISTAMINE PHOSPHATE</td>
      <td>4.0</td>
      <td>1</td>
      <td>1</td>
      <td>MOL</td>
      <td>NaN</td>
      <td>Small molecule</td>
      <td>1939.0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>NaN</td>
      <td>0.0</td>
      <td>NaN</td>
      <td>0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>111.15</td>
      <td>-0.09</td>
      <td>2.0</td>
      <td>2.0</td>
      <td>54.70</td>
      <td>2.0</td>
      <td>Y</td>
      <td>0.0</td>
      <td>NaN</td>
      <td>9.58</td>
      <td>-0.70</td>
      <td>-2.85</td>
      <td>BASE</td>
      <td>307.14</td>
      <td>1.0</td>
      <td>8.0</td>
      <td>0.56</td>
      <td>111.0796</td>
      <td>C5H15N3O8P2</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>0.0</td>
      <td>0.00</td>
    </tr>
    <tr>
      <th>4</th>
      <td>CHEMBL449</td>
      <td>CCC(C)C1(CC)C(=O)NC(=O)NC1=O</td>
      <td>2393</td>
      <td>BUTABARBITAL</td>
      <td>4.0</td>
      <td>1</td>
      <td>0</td>
      <td>MOL</td>
      <td>3228.0</td>
      <td>Small molecule</td>
      <td>1939.0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>NaN</td>
      <td>1.0</td>
      <td>-barb-</td>
      <td>0</td>
      <td>-barb-</td>
      <td>barbituric acid derivatives</td>
      <td>Sedative-Hypnotic</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>212.25</td>
      <td>0.79</td>
      <td>3.0</td>
      <td>2.0</td>
      <td>75.27</td>
      <td>3.0</td>
      <td>N</td>
      <td>0.0</td>
      <td>7.48</td>
      <td>NaN</td>
      <td>1.45</td>
      <td>1.19</td>
      <td>NEUTRAL</td>
      <td>212.25</td>
      <td>0.0</td>
      <td>15.0</td>
      <td>0.68</td>
      <td>212.1161</td>
      <td>C10H16N2O3</td>
      <td>5.0</td>
      <td>2.0</td>
      <td>0.0</td>
      <td>0.32</td>
    </tr>
  </tbody>
</table>
</div>


### **4. Read in a molecule with RDKit**

Start by reading a single SMILES in with RDKit's `MolFromSmiles()` function. This function parses the SMILES string, then constructs a graph-based molecular representation. We'll take a look at the 11th entry, sulfadiazine sodium, which is an antibacterial drug administered as a sodium salt used to treat various bacterial infections and happens to be a nice example to demonstrate the standardization process.

```python
# Read SMILES
ii = 10
smi = df["smiles"].iloc[ii]
id = df["pref_name"].iloc[ii]
mol = rdkit.Chem.MolFromSmiles(
    smi,
    sanitize = False
)
show_acs1996(mol, legend=id)
```
    
![png](rdkit_molecule_2d_profiling_files/rdkit_molecule_2d_profiling_8_0.png)
    

### **5. Process the molecule object**
Molecular standardization represents one of the most critical steps in cheminformatics. Three key processing steps are implemented here:

**Sanitization**

This step involves cleaning up the molecule and checking chemical consistency (*e.g.*, valence, kekulize, etc.). See the [RDKit book](https://www.rdkit.org/docs/RDKit_Book.html#molecular-sanitization) for details.

Without proper sanitization, downstream calculations may produce nonsensical results or fail entirely. RDKit's sanitization is particularly thorough.

**Largest Fragment Extraction** 

This step addresses an important challenge in chemical datasets: many entries represent salts or mixtures rather than single compounds. The `LargestFragmentChooser()` and `chooseInPlace()` function analyzes the molecular graph to identify disconnected components, returning the largest connected subgraph. In our example, this removes the Na ion.

**Neutralization (uncharging)**

This step standardizes the ionization state of molecules, which can vary significantly depending on the source database and pH assumptions.

This step is particularly important when calculating physicochemical descriptors that depend on ionization state, such as logP values. However, it should be applied cautiously, as bioactive conformations may require specific charge states for receptor binding.

Compare the molecule images below to the one above to see the ion and charge neutralized.

```python
# Sanitize, extract the largest fragment, and uncharge
rdkit.Chem.SanitizeMol(mol)
largest_frag_app = rdMolStandardize.LargestFragmentChooser()
uncharge_app = rdMolStandardize.Uncharger()
largest_frag_app.chooseInPlace(mol)
uncharge_app.unchargeInPlace(mol)
show_acs1996(mol, legend=id)
```
    
![png](rdkit_molecule_2d_profiling_files/rdkit_molecule_2d_profiling_10_0.png)
    

### **6. Calculate RDKit molecular descriptors**
RDKit has a lot of molecular descriptors, in fact, 217 at the writing of this post. The descriptors functions can be accessed from the `Descriptors` app in the `rdkit.Chem` module. Then looping through the `_descList` dictionary, all the descriptors can be calculated for each molecule. Subsets of descriptors can be calculated, if specific descriptors are needed.

```python
# Calculate all RDKit descriptors
descriptors = {
    "SMILES": rdkit.Chem.MolToSmiles(mol),
}
for descriptor, fxn in rdkit.Chem.Descriptors._descList:
    try:
        value = fxn(mol)
    except Exception as e:
        value = None
    descriptors[f"{descriptor}"] = value
descriptors
```

<div style="max-height: 400px; overflow-y: auto; font-family: monospace; white-space: pre;">
    {'SMILES': 'Nc1ccc(S(=O)(=O)Nc2ncccn2)cc1',
     'MaxAbsEStateIndex': np.float64(11.880950491307635),
     'MaxEStateIndex': np.float64(11.880950491307635),
     'MinAbsEStateIndex': np.float64(0.03178240740740779),
     'MinEStateIndex': np.float64(-3.6606481481481468),
     'qed': 0.7871395577774678,
     'SPS': None,
     'MolWt': 250.28300000000004,
     'HeavyAtomMolWt': 240.20299999999997,
     'ExactMolWt': 250.05244656,
     'NumValenceElectrons': 88,
     'NumRadicalElectrons': 0,
     'MaxPartialCharge': 0.2637233954888739,
     'MinPartialCharge': -0.39872771712743976,
     'MaxAbsPartialCharge': 0.39872771712743976,
     'MinAbsPartialCharge': 0.2637233954888739,
     'FpDensityMorgan1': None,
     'FpDensityMorgan2': None,
     'FpDensityMorgan3': None,
     'BCUT2D_MWHI': 32.23327104557692,
     'BCUT2D_MWLOW': 10.319170390655504,
     'BCUT2D_CHGHI': 2.1415057650483593,
     'BCUT2D_CHGLO': -2.0546834089334234,
     'BCUT2D_LOGPHI': 2.136710750473383,
     'BCUT2D_LOGPLOW': -2.1428245429362125,
     'BCUT2D_MRHI': 7.923068764373256,
     'BCUT2D_MRLOW': 0.600567336991237,
     'AvgIpc': 2.3003123655293725,
     'BalabanJ': np.float64(2.285509074399646),
     'BertzCT': 595.9628413066204,
     'Chi0': np.float64(12.303118619434356),
     'Chi0n': 8.737974215908514,
     'Chi0v': 9.55447079683624,
     'Chi1': np.float64(8.077316827624754),
     'Chi1n': 4.731902005098875,
     'Chi1v': 6.215065252693267,
     'Chi2n': 3.24693262826508,
     'Chi2v': 4.929334869674595,
     'Chi3n': 1.9818256722069898,
     'Chi3v': 3.3258965641558547,
     'Chi4n': 1.1778813725604542,
     'Chi4v': 2.135887798712802,
     'HallKierAlpha': -2.15,
     'Ipc': np.float64(7471.414563239402),
     'Kappa1': 11.338809720466916,
     'Kappa2': 4.380095906611953,
     'Kappa3': 2.6977263506132765,
     'LabuteASA': 98.57832206048985,
     'PEOE_VSA1': 5.733667477162185,
     'PEOE_VSA10': 0.0,
     'PEOE_VSA11': 0.0,
     'PEOE_VSA12': 5.948339280986494,
     'PEOE_VSA13': 10.023291153407584,
     'PEOE_VSA14': 0.0,
     'PEOE_VSA2': 0.0,
     'PEOE_VSA3': 23.10784889067544,
     'PEOE_VSA4': 0.0,
     'PEOE_VSA5': 0.0,
     'PEOE_VSA6': 0.0,
     'PEOE_VSA7': 30.33183534230805,
     'PEOE_VSA8': 18.081073417909714,
     'PEOE_VSA9': 4.895483475517775,
     'SMR_VSA1': 8.417796984328938,
     'SMR_VSA10': 21.65901670907764,
     'SMR_VSA2': 0.0,
     'SMR_VSA3': 9.967957041894415,
     'SMR_VSA4': 0.0,
     'SMR_VSA5': 4.895483475517775,
     'SMR_VSA6': 10.455762341614273,
     'SMR_VSA7': 42.725522485534206,
     'SMR_VSA8': 0.0,
     'SMR_VSA9': 0.0,
     'SlogP_VSA1': 10.455762341614273,
     'SlogP_VSA10': 11.635725555670057,
     'SlogP_VSA11': 0.0,
     'SlogP_VSA12': 0.0,
     'SlogP_VSA2': 18.385754026223353,
     'SlogP_VSA3': 10.023291153407584,
     'SlogP_VSA4': 0.0,
     'SlogP_VSA5': 0.0,
     'SlogP_VSA6': 47.62100596105198,
     'SlogP_VSA7': 0.0,
     'SlogP_VSA8': 0.0,
     'SlogP_VSA9': 0.0,
     'TPSA': None,
     'EState_VSA1': np.float64(10.023291153407584),
     'EState_VSA10': np.float64(8.417796984328938),
     'EState_VSA11': np.float64(0.0),
     'EState_VSA2': np.float64(10.84382275650427),
     'EState_VSA3': np.float64(5.687386274683562),
     'EState_VSA4': np.float64(0.0),
     'EState_VSA5': np.float64(36.6591554170726),
     'EState_VSA6': np.float64(6.06636706846161),
     'EState_VSA7': np.float64(0.0),
     'EState_VSA8': np.float64(14.690051906346504),
     'EState_VSA9': np.float64(5.733667477162185),
     'VSA_EState1': np.float64(26.025085034013607),
     'VSA_EState10': np.float64(0.0),
     'VSA_EState2': np.float64(7.672129314688839),
     'VSA_EState3': np.float64(0.0),
     'VSA_EState4': np.float64(5.982518837016492),
     'VSA_EState5': np.float64(0.03178240740740779),
     'VSA_EState6': np.float64(7.46136506530987),
     'VSA_EState7': np.float64(2.904434156378601),
     'VSA_EState8': np.float64(0.0),
     'VSA_EState9': np.float64(-3.6606481481481468),
     'FractionCSP3': 0.0,
     'HeavyAtomCount': 17,
     'NHOHCount': 3,
     'NOCount': 6,
     'NumAliphaticCarbocycles': 0,
     'NumAliphaticHeterocycles': 0,
     'NumAliphaticRings': 0,
     'NumAmideBonds': 0,
     'NumAromaticCarbocycles': 0,
     'NumAromaticHeterocycles': 0,
     'NumAromaticRings': 0,
     'NumAtomStereoCenters': 0,
     'NumBridgeheadAtoms': 0,
     'NumHAcceptors': 5,
     'NumHDonors': 2,
     'NumHeteroatoms': 7,
     'NumHeterocycles': 1,
     'NumRotatableBonds': 3,
     'NumSaturatedCarbocycles': 0,
     'NumSaturatedHeterocycles': 0,
     'NumSaturatedRings': 0,
     'NumSpiroAtoms': 0,
     'NumUnspecifiedAtomStereoCenters': 0,
     'Phi': 2.921474943674645,
     'RingCount': 2,
     'MolLogP': 0.8596000000000001,
     'MolMR': 63.69490000000001,
     'fr_Al_COO': 0,
     'fr_Al_OH': 0,
     'fr_Al_OH_noTert': 0,
     'fr_ArN': 1,
     'fr_Ar_COO': 0,
     'fr_Ar_N': 2,
     'fr_Ar_NH': 0,
     'fr_Ar_OH': 0,
     'fr_COO': 0,
     'fr_COO2': 0,
     'fr_C_O': 0,
     'fr_C_O_noCOO': 0,
     'fr_C_S': 0,
     'fr_HOCCN': 0,
     'fr_Imine': 0,
     'fr_NH0': 2,
     'fr_NH1': 1,
     'fr_NH2': 1,
     'fr_N_O': 0,
     'fr_Ndealkylation1': 0,
     'fr_Ndealkylation2': 0,
     'fr_Nhpyrrole': 0,
     'fr_SH': 0,
     'fr_aldehyde': 0,
     'fr_alkyl_carbamate': 0,
     'fr_alkyl_halide': 0,
     'fr_allylic_oxid': 0,
     'fr_amide': 0,
     'fr_amidine': 0,
     'fr_aniline': 2,
     'fr_aryl_methyl': 0,
     'fr_azide': 0,
     'fr_azo': 0,
     'fr_barbitur': 0,
     'fr_benzene': 1,
     'fr_benzodiazepine': 0,
     'fr_bicyclic': 0,
     'fr_diazo': 0,
     'fr_dihydropyridine': 0,
     'fr_epoxide': 0,
     'fr_ester': 0,
     'fr_ether': 0,
     'fr_furan': 0,
     'fr_guanido': 0,
     'fr_halogen': 0,
     'fr_hdrzine': 0,
     'fr_hdrzone': 0,
     'fr_imidazole': 0,
     'fr_imide': 0,
     'fr_isocyan': 0,
     'fr_isothiocyan': 0,
     'fr_ketone': 0,
     'fr_ketone_Topliss': 0,
     'fr_lactam': 0,
     'fr_lactone': 0,
     'fr_methoxy': 0,
     'fr_morpholine': 0,
     'fr_nitrile': 0,
     'fr_nitro': 0,
     'fr_nitro_arom': 0,
     'fr_nitro_arom_nonortho': 0,
     'fr_nitroso': 0,
     'fr_oxazole': 0,
     'fr_oxime': 0,
     'fr_para_hydroxylation': 0,
     'fr_phenol': 0,
     'fr_phenol_noOrthoHbond': 0,
     'fr_phos_acid': 0,
     'fr_phos_ester': 0,
     'fr_piperdine': 0,
     'fr_piperzine': 0,
     'fr_priamide': 0,
     'fr_prisulfonamd': 0,
     'fr_pyridine': 0,
     'fr_quatN': 0,
     'fr_sulfide': 0,
     'fr_sulfonamd': 1,
     'fr_sulfone': 0,
     'fr_term_acetylene': 0,
     'fr_tetrazole': 0,
     'fr_thiazole': 0,
     'fr_thiocyan': 0,
     'fr_thiophene': 0,
     'fr_unbrch_alkane': 0,
     'fr_urea': 0}
</div>


### **7. Loop through all molecules and merge descriptors with original dataframe**
The final two cells show how to bring the previous steps together in a loop to process all the molecules in the dataset and merge the calculated descriptors with the original dataset.  If SMILES are missing, the program skips them.

The try/except pattern is crucial because SMILES strings from public databases may fail to parse correctly. When processing thousands of molecules, robust error handling prevents pipeline failures due to individual problematic structures.

```python
# Run for all compounds
largest_frag_app = rdMolStandardize.LargestFragmentChooser()
uncharge_app = rdMolStandardize.Uncharger()
mols = []
unique_smiles = []
skipped_indices = []
for ii, smi in enumerate(df["smiles"]):
    try:
        mol = rdkit.Chem.MolFromSmiles(
            smi,
            sanitize = False
        )
        assert mol is not None
    except Exception as e:
        skipped_indices.append(ii)
        if pd.notna(smi):
            print(f"Error processing SMILES at index {ii}: {smi}")
            print(f"Exception: {e}")
        continue
    try:
        rdkit.Chem.SanitizeMol(mol)
        largest_frag_app.chooseInPlace(mol)
        uncharge_app.unchargeInPlace(mol)
    except Exception as e:
        skipped_indices.append(ii)
        print(f"Error processing molecule at index {ii}")
        print(f"Exception: {e}")
        continue
    frag_smi = rdkit.Chem.MolToSmiles(mol)
    if frag_smi in unique_smiles:
        skipped_indices.append(ii)
        continue
    descriptors = {
        "idx": ii,
        "SMILES": frag_smi,
    }
    for descriptor, fxn in Descriptors._descList:
        try:
            value = fxn(mol)
        except Exception as e:
            value = None
        descriptors[f"rdkit_{descriptor}"] = value
    mols.append(descriptors)
    unique_smiles.append(frag_smi)

print(f"{len(mols)}/{df.shape[0]} molecules processed!")
print(f"{df.shape[0]-len(mols)}/{df.shape[0]} molecules skipped!")
```
    2292/3517 molecules processed!
    1225/3517 molecules skipped!


```python
# Update dataframe with RDKit descriptors
column_update = { cc: f"chembl_{cc}" for cc in df.columns[1:] }
df = df.rename(columns=column_update)

mol_df = pd.DataFrame(mols)
mol_df.set_index("idx", inplace=True)

mol_df.columns = [ col if col == "SMILES" else col.lower() 
                   for col in mol_df.columns ]

df = pd.concat([df, mol_df], axis=1)
cols = df.columns.tolist()
cols.remove('SMILES')
cols.insert(0, 'SMILES')
df = df[cols]

# Drop skipped indices
df = df.drop(index=skipped_indices)

# Display results
print(df.shape)
df.head()
# Save to file
# df.to_csv(
#     args.output,
#     index=False,
#     encoding="utf-8"
# )
```
    (2292, 271)


<div style="overflow-x: auto; width: 100%;">
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }
    .dataframe tbody tr th {
        vertical-align: top;
    }
    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>SMILES</th>
      <th>chembl_molecule_id</th>
      <th>chembl_smiles</th>
      <th>chembl_molregno</th>
      <th>chembl_pref_name</th>
      <th>chembl_max_phase</th>
      <th>chembl_therapeutic_flag</th>
      <th>chembl_dosed_ingredient</th>
      <th>chembl_structure_type</th>
      <th>chembl_chebi_par_id</th>
      <th>chembl_molecule_type</th>
      <th>chembl_first_approval</th>
      <th>chembl_oral</th>
      <th>chembl_parenteral</th>
      <th>chembl_topical</th>
      <th>chembl_black_box_warning</th>
      <th>chembl_natural_product</th>
      <th>chembl_first_in_class</th>
      <th>chembl_chirality</th>
      <th>chembl_prodrug</th>
      <th>chembl_inorganic_flag</th>
      <th>chembl_usan_year</th>
      <th>chembl_availability_type</th>
      <th>chembl_usan_stem</th>
      <th>chembl_polymer_flag</th>
      <th>chembl_usan_substem</th>
      <th>chembl_usan_stem_definition</th>
      <th>chembl_indication_class</th>
      <th>chembl_withdrawn_flag</th>
      <th>chembl_chemical_probe</th>
      <th>chembl_orphan</th>
      <th>chembl_mw_freebase</th>
      <th>chembl_alogp</th>
      <th>chembl_hba</th>
      <th>chembl_hbd</th>
      <th>chembl_psa</th>
      <th>chembl_rtb</th>
      <th>chembl_ro3_pass</th>
      <th>chembl_num_ro5_violations</th>
      <th>chembl_cx_most_apka</th>
      <th>chembl_cx_most_bpka</th>
      <th>chembl_cx_logp</th>
      <th>chembl_cx_logd</th>
      <th>chembl_molecular_species</th>
      <th>chembl_full_mwt</th>
      <th>chembl_aromatic_rings</th>
      <th>chembl_heavy_atoms</th>
      <th>chembl_qed_weighted</th>
      <th>chembl_mw_monoisotopic</th>
      <th>chembl_full_molformula</th>
      <th>chembl_hba_lipinski</th>
      <th>chembl_hbd_lipinski</th>
      <th>chembl_num_lipinski_ro5_violations</th>
      <th>chembl_np_likeness_score</th>
      <th>rdkit_maxabsestateindex</th>
      <th>rdkit_maxestateindex</th>
      <th>rdkit_minabsestateindex</th>
      <th>rdkit_minestateindex</th>
      <th>rdkit_qed</th>
      <th>rdkit_sps</th>
      <th>rdkit_molwt</th>
      <th>rdkit_heavyatommolwt</th>
      <th>rdkit_exactmolwt</th>
      <th>rdkit_numvalenceelectrons</th>
      <th>rdkit_numradicalelectrons</th>
      <th>rdkit_maxpartialcharge</th>
      <th>rdkit_minpartialcharge</th>
      <th>rdkit_maxabspartialcharge</th>
      <th>rdkit_minabspartialcharge</th>
      <th>rdkit_fpdensitymorgan1</th>
      <th>rdkit_fpdensitymorgan2</th>
      <th>rdkit_fpdensitymorgan3</th>
      <th>rdkit_bcut2d_mwhi</th>
      <th>rdkit_bcut2d_mwlow</th>
      <th>rdkit_bcut2d_chghi</th>
      <th>rdkit_bcut2d_chglo</th>
      <th>rdkit_bcut2d_logphi</th>
      <th>rdkit_bcut2d_logplow</th>
      <th>rdkit_bcut2d_mrhi</th>
      <th>rdkit_bcut2d_mrlow</th>
      <th>rdkit_avgipc</th>
      <th>rdkit_balabanj</th>
      <th>rdkit_bertzct</th>
      <th>rdkit_chi0</th>
      <th>rdkit_chi0n</th>
      <th>rdkit_chi0v</th>
      <th>rdkit_chi1</th>
      <th>rdkit_chi1n</th>
      <th>rdkit_chi1v</th>
      <th>rdkit_chi2n</th>
      <th>rdkit_chi2v</th>
      <th>rdkit_chi3n</th>
      <th>rdkit_chi3v</th>
      <th>rdkit_chi4n</th>
      <th>rdkit_chi4v</th>
      <th>rdkit_hallkieralpha</th>
      <th>rdkit_ipc</th>
      <th>rdkit_kappa1</th>
      <th>rdkit_kappa2</th>
      <th>rdkit_kappa3</th>
      <th>rdkit_labuteasa</th>
      <th>rdkit_peoe_vsa1</th>
      <th>rdkit_peoe_vsa10</th>
      <th>rdkit_peoe_vsa11</th>
      <th>rdkit_peoe_vsa12</th>
      <th>rdkit_peoe_vsa13</th>
      <th>rdkit_peoe_vsa14</th>
      <th>rdkit_peoe_vsa2</th>
      <th>rdkit_peoe_vsa3</th>
      <th>rdkit_peoe_vsa4</th>
      <th>rdkit_peoe_vsa5</th>
      <th>rdkit_peoe_vsa6</th>
      <th>rdkit_peoe_vsa7</th>
      <th>rdkit_peoe_vsa8</th>
      <th>rdkit_peoe_vsa9</th>
      <th>rdkit_smr_vsa1</th>
      <th>rdkit_smr_vsa10</th>
      <th>rdkit_smr_vsa2</th>
      <th>rdkit_smr_vsa3</th>
      <th>rdkit_smr_vsa4</th>
      <th>rdkit_smr_vsa5</th>
      <th>rdkit_smr_vsa6</th>
      <th>rdkit_smr_vsa7</th>
      <th>rdkit_smr_vsa8</th>
      <th>rdkit_smr_vsa9</th>
      <th>rdkit_slogp_vsa1</th>
      <th>rdkit_slogp_vsa10</th>
      <th>rdkit_slogp_vsa11</th>
      <th>rdkit_slogp_vsa12</th>
      <th>rdkit_slogp_vsa2</th>
      <th>rdkit_slogp_vsa3</th>
      <th>rdkit_slogp_vsa4</th>
      <th>rdkit_slogp_vsa5</th>
      <th>rdkit_slogp_vsa6</th>
      <th>rdkit_slogp_vsa7</th>
      <th>rdkit_slogp_vsa8</th>
      <th>rdkit_slogp_vsa9</th>
      <th>rdkit_tpsa</th>
      <th>rdkit_estate_vsa1</th>
      <th>rdkit_estate_vsa10</th>
      <th>rdkit_estate_vsa11</th>
      <th>rdkit_estate_vsa2</th>
      <th>rdkit_estate_vsa3</th>
      <th>rdkit_estate_vsa4</th>
      <th>rdkit_estate_vsa5</th>
      <th>rdkit_estate_vsa6</th>
      <th>rdkit_estate_vsa7</th>
      <th>rdkit_estate_vsa8</th>
      <th>rdkit_estate_vsa9</th>
      <th>rdkit_vsa_estate1</th>
      <th>rdkit_vsa_estate10</th>
      <th>rdkit_vsa_estate2</th>
      <th>rdkit_vsa_estate3</th>
      <th>rdkit_vsa_estate4</th>
      <th>rdkit_vsa_estate5</th>
      <th>rdkit_vsa_estate6</th>
      <th>rdkit_vsa_estate7</th>
      <th>rdkit_vsa_estate8</th>
      <th>rdkit_vsa_estate9</th>
      <th>rdkit_fractioncsp3</th>
      <th>rdkit_heavyatomcount</th>
      <th>rdkit_nhohcount</th>
      <th>rdkit_nocount</th>
      <th>rdkit_numaliphaticcarbocycles</th>
      <th>rdkit_numaliphaticheterocycles</th>
      <th>rdkit_numaliphaticrings</th>
      <th>rdkit_numamidebonds</th>
      <th>rdkit_numaromaticcarbocycles</th>
      <th>rdkit_numaromaticheterocycles</th>
      <th>rdkit_numaromaticrings</th>
      <th>rdkit_numatomstereocenters</th>
      <th>rdkit_numbridgeheadatoms</th>
      <th>rdkit_numhacceptors</th>
      <th>rdkit_numhdonors</th>
      <th>rdkit_numheteroatoms</th>
      <th>rdkit_numheterocycles</th>
      <th>rdkit_numrotatablebonds</th>
      <th>rdkit_numsaturatedcarbocycles</th>
      <th>rdkit_numsaturatedheterocycles</th>
      <th>rdkit_numsaturatedrings</th>
      <th>rdkit_numspiroatoms</th>
      <th>rdkit_numunspecifiedatomstereocenters</th>
      <th>rdkit_phi</th>
      <th>rdkit_ringcount</th>
      <th>rdkit_mollogp</th>
      <th>rdkit_molmr</th>
      <th>rdkit_fr_al_coo</th>
      <th>rdkit_fr_al_oh</th>
      <th>rdkit_fr_al_oh_notert</th>
      <th>rdkit_fr_arn</th>
      <th>rdkit_fr_ar_coo</th>
      <th>rdkit_fr_ar_n</th>
      <th>rdkit_fr_ar_nh</th>
      <th>rdkit_fr_ar_oh</th>
      <th>rdkit_fr_coo</th>
      <th>rdkit_fr_coo2</th>
      <th>rdkit_fr_c_o</th>
      <th>rdkit_fr_c_o_nocoo</th>
      <th>rdkit_fr_c_s</th>
      <th>rdkit_fr_hoccn</th>
      <th>rdkit_fr_imine</th>
      <th>rdkit_fr_nh0</th>
      <th>rdkit_fr_nh1</th>
      <th>rdkit_fr_nh2</th>
      <th>rdkit_fr_n_o</th>
      <th>rdkit_fr_ndealkylation1</th>
      <th>rdkit_fr_ndealkylation2</th>
      <th>rdkit_fr_nhpyrrole</th>
      <th>rdkit_fr_sh</th>
      <th>rdkit_fr_aldehyde</th>
      <th>rdkit_fr_alkyl_carbamate</th>
      <th>rdkit_fr_alkyl_halide</th>
      <th>rdkit_fr_allylic_oxid</th>
      <th>rdkit_fr_amide</th>
      <th>rdkit_fr_amidine</th>
      <th>rdkit_fr_aniline</th>
      <th>rdkit_fr_aryl_methyl</th>
      <th>rdkit_fr_azide</th>
      <th>rdkit_fr_azo</th>
      <th>rdkit_fr_barbitur</th>
      <th>rdkit_fr_benzene</th>
      <th>rdkit_fr_benzodiazepine</th>
      <th>rdkit_fr_bicyclic</th>
      <th>rdkit_fr_diazo</th>
      <th>rdkit_fr_dihydropyridine</th>
      <th>rdkit_fr_epoxide</th>
      <th>rdkit_fr_ester</th>
      <th>rdkit_fr_ether</th>
      <th>rdkit_fr_furan</th>
      <th>rdkit_fr_guanido</th>
      <th>rdkit_fr_halogen</th>
      <th>rdkit_fr_hdrzine</th>
      <th>rdkit_fr_hdrzone</th>
      <th>rdkit_fr_imidazole</th>
      <th>rdkit_fr_imide</th>
      <th>rdkit_fr_isocyan</th>
      <th>rdkit_fr_isothiocyan</th>
      <th>rdkit_fr_ketone</th>
      <th>rdkit_fr_ketone_topliss</th>
      <th>rdkit_fr_lactam</th>
      <th>rdkit_fr_lactone</th>
      <th>rdkit_fr_methoxy</th>
      <th>rdkit_fr_morpholine</th>
      <th>rdkit_fr_nitrile</th>
      <th>rdkit_fr_nitro</th>
      <th>rdkit_fr_nitro_arom</th>
      <th>rdkit_fr_nitro_arom_nonortho</th>
      <th>rdkit_fr_nitroso</th>
      <th>rdkit_fr_oxazole</th>
      <th>rdkit_fr_oxime</th>
      <th>rdkit_fr_para_hydroxylation</th>
      <th>rdkit_fr_phenol</th>
      <th>rdkit_fr_phenol_noorthohbond</th>
      <th>rdkit_fr_phos_acid</th>
      <th>rdkit_fr_phos_ester</th>
      <th>rdkit_fr_piperdine</th>
      <th>rdkit_fr_piperzine</th>
      <th>rdkit_fr_priamide</th>
      <th>rdkit_fr_prisulfonamd</th>
      <th>rdkit_fr_pyridine</th>
      <th>rdkit_fr_quatn</th>
      <th>rdkit_fr_sulfide</th>
      <th>rdkit_fr_sulfonamd</th>
      <th>rdkit_fr_sulfone</th>
      <th>rdkit_fr_term_acetylene</th>
      <th>rdkit_fr_tetrazole</th>
      <th>rdkit_fr_thiazole</th>
      <th>rdkit_fr_thiocyan</th>
      <th>rdkit_fr_thiophene</th>
      <th>rdkit_fr_unbrch_alkane</th>
      <th>rdkit_fr_urea</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>CC(=O)OCC(=O)[C@H]1CC[C@H]2[C@@H]3CCC4=CC(=O)C...</td>
      <td>CHEMBL1200542</td>
      <td>CC(=O)OCC(=O)[C@H]1CC[C@H]2[C@@H]3CCC4=CC(=O)C...</td>
      <td>674493</td>
      <td>DESOXYCORTICOSTERONE ACETATE</td>
      <td>4.0</td>
      <td>1</td>
      <td>1</td>
      <td>MOL</td>
      <td>34671.0</td>
      <td>Small molecule</td>
      <td>1939.0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>NaN</td>
      <td>0.0</td>
      <td>-cort-; -ster-; -terone</td>
      <td>0</td>
      <td>-cort-; -ster-; -terone</td>
      <td>cortisone derivatives; steroids (androgens, an...</td>
      <td>Adrenocortical Steroid (salt-regulating)</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>372.51</td>
      <td>4.27</td>
      <td>4.0</td>
      <td>0.0</td>
      <td>60.44</td>
      <td>3.0</td>
      <td>N</td>
      <td>0.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>3.77</td>
      <td>3.77</td>
      <td>NEUTRAL</td>
      <td>372.51</td>
      <td>0.0</td>
      <td>27.0</td>
      <td>0.69</td>
      <td>372.2301</td>
      <td>C23H32O4</td>
      <td>4.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.96</td>
      <td>12.773246</td>
      <td>12.773246</td>
      <td>0.028221</td>
      <td>-0.376735</td>
      <td>0.694387</td>
      <td>43.222222</td>
      <td>372.505</td>
      <td>340.249</td>
      <td>372.230060</td>
      <td>148.0</td>
      <td>0.0</td>
      <td>0.302564</td>
      <td>-0.457860</td>
      <td>0.457860</td>
      <td>0.302564</td>
      <td>1.074074</td>
      <td>1.851852</td>
      <td>2.592593</td>
      <td>16.534689</td>
      <td>9.464840</td>
      <td>2.595559</td>
      <td>-2.576012</td>
      <td>2.661677</td>
      <td>-2.509003</td>
      <td>5.913122</td>
      <td>-0.148809</td>
      <td>2.972962</td>
      <td>1.555075</td>
      <td>707.166909</td>
      <td>19.396977</td>
      <td>16.883706</td>
      <td>16.883706</td>
      <td>12.753691</td>
      <td>10.654649</td>
      <td>10.654649</td>
      <td>9.889871</td>
      <td>9.889871</td>
      <td>8.768376</td>
      <td>8.768376</td>
      <td>7.361606</td>
      <td>7.361606</td>
      <td>-1.45</td>
      <td>1.338475e+06</td>
      <td>18.892162</td>
      <td>6.562307</td>
      <td>2.905355</td>
      <td>161.653161</td>
      <td>4.736863</td>
      <td>6.606882</td>
      <td>11.566490</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>5.969305</td>
      <td>14.383612</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>19.420579</td>
      <td>79.605471</td>
      <td>19.262465</td>
      <td>0.000000</td>
      <td>19.120475</td>
      <td>17.535795</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>34.501605</td>
      <td>72.137785</td>
      <td>6.606882</td>
      <td>11.649125</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>24.142677</td>
      <td>19.120475</td>
      <td>34.501605</td>
      <td>72.137785</td>
      <td>11.649125</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>60.44</td>
      <td>0.000000</td>
      <td>14.383612</td>
      <td>0.0</td>
      <td>35.107319</td>
      <td>29.957785</td>
      <td>44.945751</td>
      <td>12.496842</td>
      <td>0.000000</td>
      <td>6.076020</td>
      <td>13.847474</td>
      <td>4.736863</td>
      <td>5.016332</td>
      <td>0.0</td>
      <td>35.806943</td>
      <td>0.000000</td>
      <td>1.613361</td>
      <td>1.950787</td>
      <td>0.00000</td>
      <td>10.110728</td>
      <td>6.001849</td>
      <td>0.000000</td>
      <td>0.782609</td>
      <td>27.0</td>
      <td>0.0</td>
      <td>4.0</td>
      <td>4.0</td>
      <td>0.0</td>
      <td>4.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>6.0</td>
      <td>0.0</td>
      <td>4.0</td>
      <td>0.0</td>
      <td>4.0</td>
      <td>0.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>0.0</td>
      <td>3.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>4.591710</td>
      <td>4.0</td>
      <td>4.26670</td>
      <td>101.8400</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>5.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>2.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>N=C(N)N</td>
      <td>CHEMBL1200728</td>
      <td>Cl.N=C(N)N</td>
      <td>674679</td>
      <td>GUANIDINE HYDROCHLORIDE</td>
      <td>4.0</td>
      <td>1</td>
      <td>1</td>
      <td>MOL</td>
      <td>32735.0</td>
      <td>Small molecule</td>
      <td>1939.0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>NaN</td>
      <td>1.0</td>
      <td>guan-</td>
      <td>0</td>
      <td>guan-</td>
      <td>antihypertensives (guanidine derivatives)</td>
      <td>NaN</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>59.07</td>
      <td>-1.16</td>
      <td>1.0</td>
      <td>3.0</td>
      <td>75.89</td>
      <td>0.0</td>
      <td>N</td>
      <td>0.0</td>
      <td>NaN</td>
      <td>12.55</td>
      <td>-1.24</td>
      <td>-3.65</td>
      <td>BASE</td>
      <td>95.53</td>
      <td>0.0</td>
      <td>4.0</td>
      <td>0.24</td>
      <td>59.0483</td>
      <td>CH6ClN3</td>
      <td>3.0</td>
      <td>5.0</td>
      <td>0.0</td>
      <td>0.32</td>
      <td>6.055556</td>
      <td>6.055556</td>
      <td>0.333333</td>
      <td>-0.333333</td>
      <td>0.242956</td>
      <td>NaN</td>
      <td>59.072</td>
      <td>54.032</td>
      <td>59.048347</td>
      <td>24.0</td>
      <td>0.0</td>
      <td>0.182528</td>
      <td>-0.370334</td>
      <td>0.370334</td>
      <td>0.182528</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>14.880242</td>
      <td>11.139690</td>
      <td>1.520757</td>
      <td>-1.686834</td>
      <td>1.131063</td>
      <td>-2.195350</td>
      <td>5.713256</td>
      <td>1.386203</td>
      <td>0.811278</td>
      <td>2.803039</td>
      <td>26.264663</td>
      <td>3.577350</td>
      <td>2.154701</td>
      <td>2.154701</td>
      <td>1.732051</td>
      <td>0.827350</td>
      <td>0.827350</td>
      <td>0.455342</td>
      <td>0.455342</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>-0.73</td>
      <td>3.245112e+00</td>
      <td>3.270000</td>
      <td>0.710529</td>
      <td>0.173734</td>
      <td>24.104383</td>
      <td>11.467335</td>
      <td>0.000000</td>
      <td>5.959555</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>5.409284</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>5.959555</td>
      <td>5.409284</td>
      <td>0.000000</td>
      <td>11.467335</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>11.467335</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>5.959555</td>
      <td>0.000000</td>
      <td>5.409284</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>NaN</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>5.959555</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>11.467335</td>
      <td>5.409284</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>6.055556</td>
      <td>8.944444</td>
      <td>-0.333333</td>
      <td>0.00000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>4.0</td>
      <td>5.0</td>
      <td>3.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.580857</td>
      <td>0.0</td>
      <td>-1.16143</td>
      <td>16.1015</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>2.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>CCC(C)C1(CC)C(=O)NC(=O)NC1=O</td>
      <td>CHEMBL1200982</td>
      <td>CCC(C)C1(CC)C(=O)[N-]C(=O)NC1=O.[Na+]</td>
      <td>674933</td>
      <td>BUTABARBITAL SODIUM</td>
      <td>4.0</td>
      <td>1</td>
      <td>1</td>
      <td>MOL</td>
      <td>NaN</td>
      <td>Small molecule</td>
      <td>1939.0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>NaN</td>
      <td>1.0</td>
      <td>-barb-</td>
      <td>0</td>
      <td>-barb-</td>
      <td>barbituric acid derivatives</td>
      <td>Sedative-Hypnotic</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>212.25</td>
      <td>0.79</td>
      <td>3.0</td>
      <td>2.0</td>
      <td>75.27</td>
      <td>3.0</td>
      <td>N</td>
      <td>0.0</td>
      <td>7.48</td>
      <td>NaN</td>
      <td>1.45</td>
      <td>1.19</td>
      <td>NEUTRAL</td>
      <td>234.23</td>
      <td>0.0</td>
      <td>15.0</td>
      <td>0.68</td>
      <td>212.1161</td>
      <td>C10H15N2NaO3</td>
      <td>5.0</td>
      <td>2.0</td>
      <td>0.0</td>
      <td>0.32</td>
      <td>11.775694</td>
      <td>11.775694</td>
      <td>0.087407</td>
      <td>-1.094583</td>
      <td>0.679143</td>
      <td>NaN</td>
      <td>212.249</td>
      <td>196.121</td>
      <td>212.116092</td>
      <td>84.0</td>
      <td>0.0</td>
      <td>0.327632</td>
      <td>-0.276814</td>
      <td>0.327632</td>
      <td>0.276814</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>16.197656</td>
      <td>9.702763</td>
      <td>2.513510</td>
      <td>-2.368754</td>
      <td>2.352530</td>
      <td>-2.638437</td>
      <td>6.188352</td>
      <td>-0.148048</td>
      <td>2.169489</td>
      <td>3.105335</td>
      <td>291.789362</td>
      <td>11.637828</td>
      <td>9.216309</td>
      <td>9.216309</td>
      <td>6.984435</td>
      <td>5.154413</td>
      <td>5.154413</td>
      <td>3.871574</td>
      <td>3.871574</td>
      <td>3.275800</td>
      <td>3.275800</td>
      <td>2.029272</td>
      <td>2.029272</td>
      <td>-1.39</td>
      <td>1.757286e+03</td>
      <td>11.683475</td>
      <td>4.001507</td>
      <td>1.619083</td>
      <td>88.602190</td>
      <td>0.000000</td>
      <td>5.414990</td>
      <td>0.000000</td>
      <td>11.814359</td>
      <td>0.000000</td>
      <td>6.031115</td>
      <td>20.222652</td>
      <td>4.794537</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>27.192033</td>
      <td>12.338728</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>14.383612</td>
      <td>17.845474</td>
      <td>0.000000</td>
      <td>10.633577</td>
      <td>11.332897</td>
      <td>33.612855</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>10.633577</td>
      <td>4.794537</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>17.845474</td>
      <td>9.589074</td>
      <td>11.332897</td>
      <td>33.612855</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>NaN</td>
      <td>23.260464</td>
      <td>14.383612</td>
      <td>0.0</td>
      <td>5.917906</td>
      <td>12.841643</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>6.923737</td>
      <td>13.847474</td>
      <td>10.633577</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>34.486026</td>
      <td>4.313889</td>
      <td>-1.094583</td>
      <td>-1.044630</td>
      <td>-0.72465</td>
      <td>1.106968</td>
      <td>5.540314</td>
      <td>0.000000</td>
      <td>0.700000</td>
      <td>15.0</td>
      <td>2.0</td>
      <td>5.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>NaN</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>3.0</td>
      <td>2.0</td>
      <td>5.0</td>
      <td>1.0</td>
      <td>3.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>3.116767</td>
      <td>1.0</td>
      <td>0.79490</td>
      <td>53.8604</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>2.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>4.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>2.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>NCCc1c[nH]cn1</td>
      <td>CHEMBL3989520</td>
      <td>NCCc1c[nH]cn1.O=P(O)(O)O.O=P(O)(O)O</td>
      <td>2197391</td>
      <td>HISTAMINE PHOSPHATE</td>
      <td>4.0</td>
      <td>1</td>
      <td>1</td>
      <td>MOL</td>
      <td>NaN</td>
      <td>Small molecule</td>
      <td>1939.0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>NaN</td>
      <td>0.0</td>
      <td>NaN</td>
      <td>0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>111.15</td>
      <td>-0.09</td>
      <td>2.0</td>
      <td>2.0</td>
      <td>54.70</td>
      <td>2.0</td>
      <td>Y</td>
      <td>0.0</td>
      <td>NaN</td>
      <td>9.58</td>
      <td>-0.70</td>
      <td>-2.85</td>
      <td>BASE</td>
      <td>307.14</td>
      <td>1.0</td>
      <td>8.0</td>
      <td>0.56</td>
      <td>111.0796</td>
      <td>C5H15N3O8P2</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>0.0</td>
      <td>0.00</td>
      <td>5.265833</td>
      <td>5.265833</td>
      <td>0.671481</td>
      <td>0.671481</td>
      <td>0.560082</td>
      <td>9.625000</td>
      <td>111.148</td>
      <td>102.076</td>
      <td>111.079647</td>
      <td>44.0</td>
      <td>0.0</td>
      <td>0.092256</td>
      <td>-0.350904</td>
      <td>0.350904</td>
      <td>0.092256</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>14.853819</td>
      <td>10.337599</td>
      <td>1.874237</td>
      <td>-1.845221</td>
      <td>1.877082</td>
      <td>-2.017761</td>
      <td>4.933367</td>
      <td>0.934888</td>
      <td>2.185861</td>
      <td>2.500249</td>
      <td>133.393828</td>
      <td>5.819991</td>
      <td>4.593478</td>
      <td>4.593478</td>
      <td>3.931852</td>
      <td>2.609633</td>
      <td>2.609633</td>
      <td>1.599216</td>
      <td>1.599216</td>
      <td>0.955332</td>
      <td>0.955332</td>
      <td>0.575510</td>
      <td>0.575510</td>
      <td>-0.83</td>
      <td>1.071072e+02</td>
      <td>5.309470</td>
      <td>2.470712</td>
      <td>1.069114</td>
      <td>47.929986</td>
      <td>10.717646</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>4.983979</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>6.544756</td>
      <td>12.617665</td>
      <td>12.021248</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>9.967957</td>
      <td>5.733667</td>
      <td>6.420822</td>
      <td>6.544756</td>
      <td>18.218092</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>5.733667</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>16.512713</td>
      <td>6.420822</td>
      <td>0.000000</td>
      <td>5.693928</td>
      <td>12.524164</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>NaN</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>6.544756</td>
      <td>12.114750</td>
      <td>0.000000</td>
      <td>6.327320</td>
      <td>6.196844</td>
      <td>9.967957</td>
      <td>5.733667</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>6.818287</td>
      <td>0.000000</td>
      <td>6.300556</td>
      <td>0.000000</td>
      <td>0.00000</td>
      <td>4.376343</td>
      <td>0.671481</td>
      <td>0.000000</td>
      <td>0.400000</td>
      <td>8.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>2.0</td>
      <td>2.0</td>
      <td>3.0</td>
      <td>1.0</td>
      <td>2.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.639771</td>
      <td>1.0</td>
      <td>-0.08910</td>
      <td>31.3461</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>2.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>5</th>
      <td>Nc1ccc(S(=O)(=O)Nc2ccccn2)cc1</td>
      <td>CHEMBL700</td>
      <td>Nc1ccc(S(=O)(=O)Nc2ccccn2)cc1</td>
      <td>32842</td>
      <td>SULFAPYRIDINE</td>
      <td>4.0</td>
      <td>1</td>
      <td>1</td>
      <td>MOL</td>
      <td>132842.0</td>
      <td>Small molecule</td>
      <td>1939.0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>NaN</td>
      <td>0.0</td>
      <td>sulfa-</td>
      <td>0</td>
      <td>sulfa-</td>
      <td>antimicrobials (sulfonamides derivatives)</td>
      <td>Suppressant (dermatitis herpetiformis)</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>249.30</td>
      <td>1.46</td>
      <td>4.0</td>
      <td>2.0</td>
      <td>85.08</td>
      <td>3.0</td>
      <td>N</td>
      <td>0.0</td>
      <td>6.24</td>
      <td>2.14</td>
      <td>1.01</td>
      <td>0.24</td>
      <td>ACID</td>
      <td>249.30</td>
      <td>2.0</td>
      <td>17.0</td>
      <td>0.81</td>
      <td>249.0572</td>
      <td>C11H11N3O2S</td>
      <td>5.0</td>
      <td>3.0</td>
      <td>0.0</td>
      <td>-1.63</td>
      <td>11.920950</td>
      <td>11.920950</td>
      <td>0.150327</td>
      <td>-3.598148</td>
      <td>0.806426</td>
      <td>11.058824</td>
      <td>249.295</td>
      <td>238.207</td>
      <td>249.057198</td>
      <td>88.0</td>
      <td>0.0</td>
      <td>0.262547</td>
      <td>-0.398728</td>
      <td>0.398728</td>
      <td>0.262547</td>
      <td>1.058824</td>
      <td>1.764706</td>
      <td>2.352941</td>
      <td>32.233271</td>
      <td>10.319152</td>
      <td>2.133279</td>
      <td>-2.063880</td>
      <td>2.146726</td>
      <td>-2.131944</td>
      <td>7.923180</td>
      <td>0.600815</td>
      <td>2.300312</td>
      <td>2.285509</td>
      <td>594.373487</td>
      <td>12.303119</td>
      <td>8.868111</td>
      <td>9.684607</td>
      <td>8.077317</td>
      <td>4.872105</td>
      <td>6.355268</td>
      <td>3.389512</td>
      <td>5.071915</td>
      <td>2.094227</td>
      <td>3.464862</td>
      <td>1.266702</td>
      <td>2.275017</td>
      <td>-2.08</td>
      <td>7.471415e+03</td>
      <td>11.406730</td>
      <td>4.423184</td>
      <td>2.729239</td>
      <td>99.358674</td>
      <td>5.733667</td>
      <td>5.817863</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>10.023291</td>
      <td>0.000000</td>
      <td>4.722095</td>
      <td>13.401776</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>6.066367</td>
      <td>36.398202</td>
      <td>11.884230</td>
      <td>4.895483</td>
      <td>8.417797</td>
      <td>21.528540</td>
      <td>0.000000</td>
      <td>4.983979</td>
      <td>0.000000</td>
      <td>4.895483</td>
      <td>10.455762</td>
      <td>48.661413</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>10.455762</td>
      <td>11.505249</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>13.401776</td>
      <td>10.023291</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>53.556897</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>85.08</td>
      <td>10.023291</td>
      <td>8.417797</td>
      <td>0.0</td>
      <td>10.713346</td>
      <td>5.687386</td>
      <td>0.000000</td>
      <td>30.462312</td>
      <td>18.199101</td>
      <td>0.000000</td>
      <td>9.706073</td>
      <td>5.733667</td>
      <td>26.216196</td>
      <td>0.0</td>
      <td>4.042339</td>
      <td>0.000000</td>
      <td>6.010490</td>
      <td>0.281782</td>
      <td>10.94929</td>
      <td>1.514717</td>
      <td>0.000000</td>
      <td>-3.598148</td>
      <td>0.000000</td>
      <td>17.0</td>
      <td>3.0</td>
      <td>5.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>2.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>4.0</td>
      <td>2.0</td>
      <td>6.0</td>
      <td>1.0</td>
      <td>3.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>2.967886</td>
      <td>2.0</td>
      <td>1.46460</td>
      <td>65.8999</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>2.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
  </tbody>
</table>
</div>


## Conclusion

In this tutorial, we've covered essential techniques for standardizing and analyzing chemical structures with RDKit. From SMILES processing to molecular standardization (sanitization, salt removal, and charge neutralization) to descriptor calculation. These approaches ensure consistency across chemical datasets, enabling reliable property predictions and comparisons critical for drug discovery applications.
