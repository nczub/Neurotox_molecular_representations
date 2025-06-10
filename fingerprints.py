import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys, rdMolDescriptors


def calculate_fingerprints(smiles):
    fingerprints = []
    for smi in smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            # Morgan Fingerprint (radius 2, 2048 bits)
            morgan_fp = list(AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048))

            # MACCS Fingerprint (166 bits)
            maccs_fp = list(MACCSkeys.GenMACCSKeys(mol))

            # RDKit Fingerprint (2048 bits)
            rdkit_fp = list(rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol, nBits=2048))

            fingerprints.append((morgan_fp, maccs_fp, rdkit_fp))
        else:
            fingerprints.append((None, None, None))
    return fingerprints


def main(input_csv, output_morgan_csv, output_maccs_csv, output_rdkit_csv):
    data = pd.read_csv(input_csv)

    if 'SMILES' not in data.columns:
        raise ValueError("CSV file must have column 'SMILES'.")

    smiles = data['SMILES']

    
    fingerprints = calculate_fingerprints(smiles)

    
    morgan_results = []
    maccs_results = []
    rdkit_results = []

    for idx, (morgan, maccs, rdkit) in enumerate(fingerprints):
        smi = smiles[idx]
        if morgan is not None:
            morgan_results.append({'SMILES': smi, **{f'Morgan_FP_{i}': bit for i, bit in enumerate(morgan)}})
            maccs_results.append({'SMILES': smi, **{f'MACCS_FP_{i}': bit for i, bit in enumerate(maccs)}})
            rdkit_results.append({'SMILES': smi, **{f'RDKit_FP_{i}': bit for i, bit in enumerate(rdkit)}})
        else:
            morgan_results.append({'SMILES': smi})
            maccs_results.append({'SMILES': smi})
            rdkit_results.append({'SMILES': smi})


    pd.DataFrame(morgan_results).to_csv(output_morgan_csv, index=False)
    pd.DataFrame(maccs_results).to_csv(output_maccs_csv, index=False)
    pd.DataFrame(rdkit_results).to_csv(output_rdkit_csv, index=False)

    print(f"Morgan fingerprins saved in file: {output_morgan_csv}")
    print(f"MACCS fingerprins saved in file: {output_maccs_csv}")
    print(f"RDKit fingerprins saved in file: {output_rdkit_csv}")

if __name__ == "__main__":
    input_file = f"neurotox_MEA_SMILES_part.csv"
    output_morgan_file = "neurotox_MEA_with_SMILES_morgan_fingerprints.csv"
    output_maccs_file = "neurotox_MEA_with_SMILES_maccs_fingerprints.csv"
    output_rdkit_file = "neurotox_MEA_with_SMILES_rdkit_fingerprints.csv"

    main(input_file, output_morgan_file, output_maccs_file, output_rdkit_file)



# add parameters of experiments and delete SMILES column
df_maccs = pd.read_csv("neurotox_MEA_with_SMILES_maccs_fingerprints.csv")
df_morgan = pd.read_csv("neurotox_MEA_with_SMILES_morgan_fingerprints.csv")
df_rdkit = pd.read_csv("neurotox_MEA_with_SMILES_rdkit_fingerprints.csv")

df_core = pd.read_csv(f"neurotox_MEA_SMILES_part.csv")

data_all_morgan = pd.concat([df_core.iloc[:, :9], df_morgan, df_core.iloc[:, -1]], axis = 1)
del data_all_morgan['SMILES']
data_all_morgan.to_csv("neurotox_MEA_morgan_fingerprints.csv", index = False, header = True)

data_all_maccs = pd.concat([df_core.iloc[:, :9], df_maccs, df_core.iloc[:, -1]], axis = 1)
del data_all_maccs['SMILES']
data_all_maccs.to_csv("neurotox_MEA_maccs_fingerprints.csv", index = False, header = True)

data_all_rdkit = pd.concat([df_core.iloc[:, :9], df_rdkit, df_core.iloc[:, -1]], axis = 1)
del data_all_rdkit['SMILES']
data_all_rdkit.to_csv("neurotox_MEA_rdkit_fingerprints.csv", index = False, header = True)