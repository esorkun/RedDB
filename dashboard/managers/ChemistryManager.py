# -*- coding: utf-8 -*-
from rdkit.Chem import AllChem as rdkit


class ChemistryManager():
    
    def GenerateInchiKeyFromSmiles(smiles):
        return rdkit.inchi.MolToInchiKey(ChemistryManager.GenerateMolFromSmiles(smiles))
    
    def GenerateMolFromSmiles(smiles):
        return rdkit.MolFromSmiles(smiles)
        
