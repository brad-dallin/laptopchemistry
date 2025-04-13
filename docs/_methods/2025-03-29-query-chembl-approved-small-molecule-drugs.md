---
layout: post
title: "How to Query ChEMBL Database for Approved Small Molecule Drugs: A Step-by-Step Guide"
author: Brad Dallin
catagories: methods
---

In this step-by-step guide, I'll show you how to build a SQL query to extract chemical structures and drug data from the ChEMBL database. As a valuable resource for drug discovery and cheminformatics research, ChEMBL contains extensive information on bioactive molecules. This tutorial assumes some basic familiarity with SQL and Python.

## **0. Setting Up Your Environment**

Before querying the ChEMBL database, you'll need to set up your environment properly. This section covers the prerequisite installations and configurations needed to follow this tutorial.


### **Prerequisites**

ChEMBL Database: A local installation of the ChEMBL database

PostgreSQL: A running PostgreSQL server (I used Postgres App for Mac)

Python 3.7+: With necessary libraries for database connection and data manipulation


### **Installing the ChEMBL Database**

The ChEMBL database is a comprehensive resource of bioactive molecules with drug-like properties. To download and install:

Install the Postgres App ([here](https://postgresapp.com)), if you're using macOS

```python
# Install PostgreSQL (v17 - 2025-03-01)
curl -JOfSL --output-dir {DOWNLOAD_PATH} https://github.com/PostgresApp/PostgresApp/releases/download/v2.8.1/Postgres-2.8.1-17.dmg
open Postgres-2.8.1-17.dmg

# Add CLI utilities
sudo mkdir -p /etc/paths.d && echo /Applications/Postgres.app/Contents/Versions/latest/bin | sudo tee /etc/paths.d/postgresapp
```

Download the latest PostgreSQL dump from the [ChEMBL Downloads page](https://chembl.gitbook.io/chembl-interface-documentation/downloads)

```python
# Download Chembl DB (v35 - 2025-03-01)
curl -JOfSL --output-dir {DOWNLOAD_PATH} https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_35_postgresql.tar.gz
tar -xzvf chembl_35_postgresql.tar.gz
```

### **Create a new database for ChEMBL**
This process can take some time depending on your system specifications.

```python
# Create new database
dropdb chembl35
createdb chembl35
pg_restore --verbose --no-owner --host=localhost --port=5432 --username=postgres --dbname=chembl35 {PATHTO}/chembl_35_postgresql.dmp
```


### **Python Environment Setup**

For this tutorial, we'll use the adbc_driver_postgresql module for database connectivity and pandas for data manipulation. Install these packages using conda, preferably in a clean environment. I installed these within my RDKit environment.

```python
# Create RDKit conda environment
conda create -n rdkit-env conda-forge::rdkit

#Conda install ADBC Driver PostgreSQL
conda install conda-forge::adbc-driver-postgresql

# Conda install PyArrow
conda install conda-forge::pyarrow
```

The psycopg2 module can also be used, if you prefer. Although, the adbc_driver_postgresql module provides a more modern interface compared to psycopg2 and integrates well with pandas for data analysis workflow, see [here](https://pandas.pydata.org/docs/reference/api/pandas.read_sql.html).

Now that our environment is ready, let's proceed to connect to the database and start exploring the ChEMBL data structure.


## **1. Import modules**

Import pandas and adbc_driver_postgresql. By default, pandas truncates the number of columns printed. I adjusted the max column setting to print all columns.

```python
# Import modules
import pandas as pd
import adbc_driver_postgresql
import adbc_driver_postgresql.dbapi

# Expand to see all columns
pd.set_option('display.max_columns', None)

# Print versions
print(f"Pandas Version: {pd.__version__}")
print(f"ADBC Driver Version: {adbc_driver_postgresql.__version__}")
```
    Pandas Version: 2.2.3
    ADBC Driver Version: 1.5.0


## 2. **Test Connection to PostgreSQL**

Set the URI hosting the ChEMBL database to localhost. Then use a basic SQL query to test the connection to the database.

```python
uri = "postgresql://localhost/chembl35"
try:
    conn = adbc_driver_postgresql.dbapi.connect(uri)
    with conn.cursor() as cur:
        cur.execute("SELECT 1")
        assert cur.fetchone() == (1,)
    print("Database connected successfully")
except:
    print("Database not connected successfully")
```
    Database connected successfully


## 3. **Preview All Available Tables in ChEMBL Database**

The ChEMBL database contains a lot of data which has been split into tables. You will often need to query and merge multiple tables to get all the data needed. Run the query below to list all the available tables in the database.

```python
sql_list_tables = """
SELECT table_name
FROM information_schema.tables
    WHERE table_schema = 'public'
ORDER BY table_name
"""
with conn.cursor() as cur:
    cur.execute(sql_list_tables)
    chembl_tables_df = pd.DataFrame(
        cur.fetchall(),
        columns=[desc[0] for desc in cur.description]
    )
print(chembl_tables_df["table_name"].values)
```
    ['action_type' 'activities' 'activity_properties' 'activity_smid'
     'activity_stds_lookup' 'activity_supp' 'activity_supp_map'
     'assay_class_map' 'assay_classification' 'assay_parameters' 'assay_type'
     'assays' 'atc_classification' 'binding_sites' 'bio_component_sequences'
     'bioassay_ontology' 'biotherapeutic_components' 'biotherapeutics'
     'cell_dictionary' 'chembl_id_lookup' 'chembl_release' 'component_class'
     'component_domains' 'component_go' 'component_sequences'
     'component_synonyms' 'compound_properties' 'compound_records'
     'compound_structural_alerts' 'compound_structures'
     'confidence_score_lookup' 'curation_lookup' 'data_validity_lookup'
     'defined_daily_dose' 'docs' 'domains' 'drug_indication' 'drug_mechanism'
     'drug_warning' 'formulations' 'frac_classification' 'go_classification'
     'hrac_classification' 'indication_refs' 'irac_classification'
     'ligand_eff' 'mechanism_refs' 'metabolism' 'metabolism_refs'
     'molecule_atc_classification' 'molecule_dictionary'
     'molecule_frac_classification' 'molecule_hierarchy'
     'molecule_hrac_classification' 'molecule_irac_classification'
     'molecule_synonyms' 'organism_class' 'patent_use_codes'
     'predicted_binding_domains' 'product_patents' 'products'
     'protein_class_synonyms' 'protein_classification' 'relationship_type'
     'research_companies' 'research_stem' 'site_components' 'source'
     'structural_alert_sets' 'structural_alerts' 'target_components'
     'target_dictionary' 'target_relations' 'target_type' 'tissue_dictionary'
     'usan_stems' 'variant_sequences' 'version' 'warning_refs']


## 4. **Preview All Columns in Compound Structures Table**

Run this query to see which columns are available in the compound structures table. Canonical SMILES are also available in this table.

```python
sql_table = "compound_structures"
sql_list_columns = """
SELECT column_name
FROM information_schema.columns
    WHERE table_name = 'compound_records'
ORDER BY ordinal_position
"""
with conn.cursor() as cur:
    cur.execute(sql_list_columns)
    chembl_columns_df = pd.DataFrame(
        cur.fetchall(),
        columns=[desc[0] for desc in cur.description]
    )
print(f"Table name: {sql_table}")
print(f"Columns: {chembl_columns_df["column_name"].values}")
```
    Table name: compound_structures
    Columns: ['record_id' 'molregno' 'doc_id' 'compound_key' 'compound_name' 'src_id'
     'src_compound_id' 'cidx']


## 5. **Preview All Columns in Molecule Dictionary**

Run this query to see which columns are available in the molecule dictionary table, I ran this query. The molecule dictionary table contains much of the drug information data.

```python
sql_table = "molecule_dictionary"
sql_list_columns = """
SELECT column_name
FROM information_schema.columns
    WHERE table_name = 'molecule_dictionary'
ORDER BY ordinal_position
"""
with conn.cursor() as cur:
    cur.execute(sql_list_columns)
    chembl_columns_df = pd.DataFrame(
        cur.fetchall(),
        columns=[desc[0] for desc in cur.description]
    )
print(f"Table name: {sql_table}")
print(f"Columns: {chembl_columns_df["column_name"].values}")
```
    Table name: molecule_dictionary
    Columns: ['molregno' 'pref_name' 'chembl_id' 'max_phase' 'therapeutic_flag'
     'dosed_ingredient' 'structure_type' 'chebi_par_id' 'molecule_type'
     'first_approval' 'oral' 'parenteral' 'topical' 'black_box_warning'
     'natural_product' 'first_in_class' 'chirality' 'prodrug' 'inorganic_flag'
     'usan_year' 'availability_type' 'usan_stem' 'polymer_flag' 'usan_substem'
     'usan_stem_definition' 'indication_class' 'withdrawn_flag'
     'chemical_probe' 'orphan']


## 6. **Preview All Columns in Compound Properties**

Run this query to see which columns are available in the compound properties table. Compound properties are mostly the physicochemical properties of the molecule.

```python
sql_table = "compound_properties"
sql_list_columns = """
SELECT column_name
FROM information_schema.columns
    WHERE table_name = 'compound_properties'
ORDER BY ordinal_position
"""
with conn.cursor() as cur:
    cur.execute(sql_list_columns)
    chembl_columns_df = pd.DataFrame(
        cur.fetchall(),
        columns=[desc[0] for desc in cur.description]
    )
print(f"Table name: {sql_table}")
print(f"Columns: {chembl_columns_df["column_name"].values}")
```
    Table name: compound_properties
    Columns: ['molregno' 'mw_freebase' 'alogp' 'hba' 'hbd' 'psa' 'rtb' 'ro3_pass'
     'num_ro5_violations' 'cx_most_apka' 'cx_most_bpka' 'cx_logp' 'cx_logd'
     'molecular_species' 'full_mwt' 'aromatic_rings' 'heavy_atoms'
     'qed_weighted' 'mw_monoisotopic' 'full_molformula' 'hba_lipinski'
     'hbd_lipinski' 'num_lipinski_ro5_violations' 'np_likeness_score']


## 7. **Query ChEMBL Database for Small Molecule Drugs and Merge Tables**

You can take the information about the columns above and combine it into a query/filter for small molecule drugs. It may be a little excessive to include all columns, but this is how you could keep all the data in view. It’s often easier to filter data later, rather than go back to retrieve or generate the data again. Uncomment chembl_df.to_csv command and add a file path to save the data to file.

```python
sql_query="""
SELECT DISTINCT m.chembl_id AS molecule_chembl_id,
s.canonical_smiles AS smiles,
-- Molecule Dictionary
m.molregno, m.pref_name, m.max_phase, m.therapeutic_flag, m.dosed_ingredient,
m.structure_type, m.chebi_par_id, m.molecule_type, m.first_approval, m.oral,
m.parenteral, m.topical, m.black_box_warning, m.natural_product, m.first_in_class,
m.chirality, m.prodrug, m.inorganic_flag, m.usan_year, m.availability_type,
m.usan_stem, m.polymer_flag, m.usan_substem, m.usan_stem_definition,
m.indication_class, m.withdrawn_flag, m.chemical_probe, m.orphan,
-- Compound Properties
p.mw_freebase, p.alogp, p.hba, p.hbd, p.psa, p.rtb, p.ro3_pass,
p.num_ro5_violations, p.cx_most_apka, p.cx_most_bpka, p.cx_logp,
p.cx_logd, p.molecular_species, p.full_mwt, p.aromatic_rings,
p.heavy_atoms, p.qed_weighted, p.mw_monoisotopic, p.full_molformula,
p.hba_lipinski, p.hbd_lipinski, p.num_lipinski_ro5_violations,
p.np_likeness_score
-- Join tables
FROM compound_structures s
    RIGHT JOIN compound_properties p ON s.molregno = p.molregno
    RIGHT JOIN molecule_dictionary m ON p.molregno = m.molregno
    JOIN compound_records r ON m.molregno = r.molregno
    AND m.max_phase = 4
    AND m.molecule_type = 'Small molecule'
    AND m.inorganic_flag = 0
ORDER BY m.first_approval ASC
"""
with conn.cursor() as cur:
    cur.execute(sql_query)
    chembl_df = pd.DataFrame(
        cur.fetchall(),
        columns=[desc[0] for desc in cur.description]
    )
print(chembl_df.shape)
chembl_df.head()
# Save to file
# chembl_df.to_csv(
#     output_file_path,
#     index=False,
#     encoding="utf-8"
# )
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
      <th>molecule_chembl_id</th>
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
      <td>None</td>
      <td>None</td>
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
      <td>None</td>
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
      <td>None</td>
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
      <td>None</td>
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
      <td>None</td>
      <td>0</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
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
      <td>None</td>
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
      <td>None</td>
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
