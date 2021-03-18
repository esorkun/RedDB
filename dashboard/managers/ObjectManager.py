# -*- coding: utf-8 -*-


from .MoleculeManager import MoleculeManager

class FileFormat(): 
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
            
class PairSmilesFileFormat(): 
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
            
class PairInfo():
    def FormatFileName(name):
        actualName=""
        if ("frag" in name):
            actualName = name.replace(";", "").replace("[", "").replace("]", "").replace(" ", "")
        elif ("_-_" in name):
            actualName = name.replace(",", "").replace(" ", "").split("[")[0]
        else:
            print("ERROR")
        return actualName
    
    def GenerateMolFromSmiles(smiles):
        return MoleculeManager.MolFromSmiles(smiles)
    
    def GenerateInchiKeyFromSmiles(smiles):
        return MoleculeManager.InchiKeyFromSmiles(smiles)
    
    ct = 0
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
            if (key == "name"):
                setattr(self, "actualName", PairInfo.FormatFileName(value))
            elif (key == "smiles"):
                mol = PairInfo.GenerateMolFromSmiles(value)
                setattr(self, "isParent", None)
                setattr(self, "mol", mol)
                setattr(self, "inchiKey", PairInfo.GenerateInchiKeyFromSmiles(value))
       
class MoleculeInfo():
    def FormatFileName(name):
        actualName=""
        if ("frag" in name):
            actualName = name.replace(";", "").replace("[", "").replace("]", "").replace(" ", "")
        elif ("_-_" in name):
            actualName = name.replace(",", "").replace(" ", "").split("[")[0]
        else:
            print("ERROR")
        return actualName
    
    def GenerateMolFromSmiles(smiles):
        return MoleculeManager.MolFromSmiles(smiles)
    
    def GenerateInchiKeyFromSmiles(smiles):
        return MoleculeManager.InchiKeyFromSmiles(smiles)
    
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
            if (key == "name"):
                setattr(self, "actualName", MoleculeInfo.FormatFileName(value))
            elif (key == "smiles"):
                mol = MoleculeInfo.GenerateMolFromSmiles(value)
                setattr(self, "isParent", None)
                setattr(self, "mol", mol)
                setattr(self, "inchiKey", MoleculeInfo.GenerateInchiKeyFromSmiles(value))
                
class JobInfo():
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        
        