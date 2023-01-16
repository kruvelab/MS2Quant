
from rdkit import Chem
from rdkit.Chem import AllChem

import pandas


#nBits = 1024
#radius = 2
#m1 = Chem.MolFromSmiles('N#CC(OC1OC(COC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O)c1ccccc1')
#fp1 = AllChem.GetMorganFingerprintAsBitVect(m1, radius = radius, nBits = nBits)
#fp1_name = [f'Bit_{i}' for i in range(nBits)]
#fp1_bits = [list(l) for l in fp1]
#print(fp1_bits)
#<rdkit.DataStructs.cDataStructs.UIntSparseIntVect object at 0x...>
#m2 = Chem.MolFromSmiles('Cc1ncccc1')
#fp2 = AllChem.GetMorganFingerprint(m2,2)
#DataStructs.DiceSimilarity(fp1,fp2)
#0.55...

#fp1 = AllChem.GetMorganFingerprintAsBitVect(m1,2,nBits=1024)
#fp1
#<rdkit.DataStructs.cDataStructs.ExplicitBitVect object at 0x...>
#fp2 = AllChem.GetMorganFingerprintAsBitVect(m2,2,nBits=1024)
#DataStructs.DiceSimilarity(fp1,fp2)
#0.51...

def Morgan2_calculation(SMILES):
    mol = Chem.MolFromSmiles(SMILES)
    #info = {}
    SMILES_descriptors = AllChem.GetMorganFingerprintAsBitVect(mol, 
                                                               radius = 2, 
                                                               nBits = 1024)
    #SMILES_descriptors = SMILES_descriptors.ToBinary()
    SMILES_descriptors = list(SMILES_descriptors)
    SMILES_descriptors = pandas.DataFrame(SMILES_descriptors)

    return SMILES_descriptors


test = Morgan2_calculation('N#CC(OC1OC(COC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O)c1ccccc1')
#print(test)

#ecfp6_name = [f'Bit_{i}' for i in range(nBits)]
#ecfp6_bits = [list(l) for l in ECFP6]
#df_morgan = pd.DataFrame(ecfp6_bits, index = esol_data.smiles, columns=ecfp6_name)
#df_morgan.head(1)
#test.ToBinary()
#print(test)