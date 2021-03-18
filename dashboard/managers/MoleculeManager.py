# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 00:42:56 2020

@author: Elif
"""

from rdkit.Chem import AllChem
import rdkit as rdkit

class MoleculeManager():
    
    def InchiKeyFromSmiles(smiles):
        mol = AllChem.MolFromSmiles(smiles)
        InchiKey = AllChem.inchi.MolToInchiKey(mol)
        return InchiKey
    
    def InchiKeyFromMol(mol):
        InchiKey = AllChem.inchi.MolToInchiKey(mol)
        return InchiKey
    
    def MolFromSmiles(smiles):
        return rdkit.Chem.MolFromSmiles(smiles)
