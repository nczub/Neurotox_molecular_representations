import pandas as pd
from rdkit import Chem
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import rdPartialCharges

def calculate_descriptors(smiles_col):
    descriptor_names = [desc[0] for desc in Descriptors._descList]
    calculator = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names)
    descriptors_data = []

    for smiles in smiles_col:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            descriptors = calculator.CalcDescriptors(mol)
            descriptors_data.append(descriptors)
        else:
            descriptors_data.append([None] * len(descriptor_names))

    descriptors_df = pd.DataFrame(descriptors_data, columns=descriptor_names)
    return descriptors_df

file_name ="neurotox_MEA_SMILES_part.csv"
df = pd.read_csv(file_name)
parameters = df.iloc[:, :9]
df_output = df.iloc[:, -1]
descriptors_df = calculate_descriptors(df["SMILES"])

result_df = pd.concat([parameters, descriptors_df, df_output], axis=1)
del result_df['SMILES']
result_df.to_csv(f"neurotox_MEA_rdkit_descriptors.csv", index=False)