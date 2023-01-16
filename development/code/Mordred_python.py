# -*- coding: utf-8 -*-
#Data preprocessing

# Import the libraries


from rdkit import Chem
from mordred import Calculator, descriptors
import pandas


def Mordred_calculation(SMILES):
    calc = Calculator(descriptors, ignore_3D=True)
    mol = Chem.MolFromSmiles(SMILES)
    SMILES_descriptors = calc(mol)
    return SMILES_descriptors

