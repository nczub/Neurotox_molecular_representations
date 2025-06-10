import pandas as pd
import numpy as np
from rdkit import Chem
from mordred import Calculator, descriptors

# input data with SMILES column
input_file = "neurotox_MEA_SMILES_part.csv"
input_data = pd.read_csv(input_file)

smiles_list = input_data["SMILES"]
output_variable = input_data["Mean Spike Rate (mean firing rate) [% of control]"]
data_id = input_data.iloc[:, :-1]

calc = Calculator(descriptors, ignore_3D=True)

descriptors = []
valid_smiles = []
valid_indices = []

for i, smi in enumerate(smiles_list):
    mol = Chem.MolFromSmiles(smi)
    if mol is not None:
        valid_smiles.append(smi)
        valid_indices.append(i)
        desc_df = calc.pandas([mol])
        descriptors.append(desc_df)

descriptor_df = pd.concat(descriptors, ignore_index=True)
data_id_valid = data_id.iloc[valid_indices].reset_index(drop=True)
output_variable_valid = output_variable.iloc[valid_indices].reset_index(drop=True)

combined_df = pd.concat([data_id_valid, descriptor_df, output_variable_valid], axis=1)
temp_file = "neurotox_MEA_mordred.csv"
combined_df.to_csv(temp_file, index=False)

# Data cleaning

df = pd.read_csv(temp_file)

parameters = df.iloc[:, :9]
descriptors_clean = df.iloc[:, 9:-1]

# change error values into NaN
descriptors_clean = descriptors_clean.applymap(lambda x: x if not isinstance(x, str) else np.nan)

# fill empty cells with mean value
descriptors_clean.fillna(descriptors_clean.mean(), inplace=True)

for col in ["Lipinski", "GhoseFilter"]:
    if col in descriptors_clean.columns:
        descriptors_clean[col] = descriptors_clean[col].astype(int)

descriptors_clean.fillna(0, inplace=True)
# remove columns with constant values
descriptors_clean = descriptors_clean.loc[:, (descriptors_clean != 0).any(axis=0)]

cleaned_data = pd.concat([parameters, descriptors_clean, df.iloc[:, -1]], axis=1)
del cleaned_data['SMILES']

output_file = "neurotox_MEA_mordred_clean.csv"
cleaned_data.to_csv(output_file, index=False)