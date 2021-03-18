# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 21:43:05 2021

@author: Elif
"""




import math
from rdkit import Chem
from rdkit.Chem import rdMolHash as rdMolHash
import rdkit
from rdkit.Chem import rdMolDescriptors

from rdkit.Chem.rdMolDescriptors import CalcMolFormula       


import pandas as p
import os
path = os.path.join(os.getcwd(), 'C:\\Users\\Elif\\Desktop\\RedDB_molecule.csv')
df = p.read_csv(path)
print("DF len : " + str(len(df)))
count=0

def isnan(value):
    try:
        return math.isnan(float(value))
    except:
        return False

def GetNumber(stoichiometry, atom):
    num=""
    splitedFormula =""
    splitedFormula = stoichiometry.split(atom)
    if len(splitedFormula)== 1:
        num=1
    else :
        secondPart = splitedFormula[1]
        checkIsItFirst = False
        for val in secondPart:
            if val.isalpha():
                break;
            else: 
                checkIsItFirst = True
                num += str(val)
                secondPart=secondPart[1:]
        if not checkIsItFirst:
            num=1
    return num

for index, row in df.iterrows():
    smiles = row['smiles']
    readedFormula = row['stoichiometry']
    molObj = Chem.MolFromSmiles(smiles)
    formula = CalcMolFormula(molObj)
    atoms = rdkit.Chem.rdchem.Mol.GetAtoms(molObj)
    if isnan(formula) or isnan(readedFormula):
        print("*************** NAN VALUE : " + str(readedFormula) + " - " + str(formula))
    else:
        if readedFormula != formula:
            for atom in atoms:
                atomSymbol = rdkit.Chem.rdchem.Atom.GetSymbol(atom)
                atomSymbolSTR = str(atomSymbol)
                readedNumber = GetNumber(str(readedFormula), atomSymbolSTR)
                createdNumber = GetNumber(str(formula), atomSymbolSTR)
                if readedNumber != createdNumber:
                    print(readedFormula)
                    print(formula)
                    print(count)
                    print("-----------------------")
                    count+=1
                    break
    
# print("count : "+ str(count))
