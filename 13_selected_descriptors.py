import pandas as pd
from rdkit import Chem
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import rdPartialCharges

def calculate_selected_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    
    descriptors = {
        "MolWt": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "MolMR": Descriptors.MolMR(mol),
        "TPSA": Descriptors.TPSA(mol),
        "NumHDonors": Descriptors.NumHDonors(mol),
        "NumHAcceptors": Descriptors.NumHAcceptors(mol),
        "NumAromaticRings": Descriptors.NumAromaticRings(mol),
        "NumAliphaticRings": Descriptors.NumAliphaticRings(mol),
        "NumRotatableBonds": Descriptors.NumRotatableBonds(mol),
        "FractionCSP3": Descriptors.FractionCSP3(mol),
        "NumRadicalElectrons": Descriptors.NumRadicalElectrons(mol),
        "NumSaturatedRings": Descriptors.NumSaturatedRings(mol),
    }
    nitro_smarts = Chem.MolFromSmarts("[NX3](=O)[O-]")
    descriptors["NumNitroGroups"] = len(mol.GetSubstructMatches(nitro_smarts))
    return descriptors


file_name ="neurotox_MEA_SMILES_part.csv"
df = pd.read_csv(file_name)
parameters = df.iloc[:, :9]
df_output = df.iloc[:, -1]
results = [calculate_selected_descriptors(smiles) for smiles in df["SMILES"]]
descriptors_df = pd.DataFrame(results)

final_df = pd.concat([parameters, descriptors_df, df_output], axis=1)
del final_df['SMILES']
final_df.to_csv(f"neurotox_MEA_13_descriptors.csv", index=False, header = True)