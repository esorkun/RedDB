

from rdkit import Chem
from rdkit.Chem import rdMolHash as rdMolHash
import rdkit
from rdkit.Chem import rdMolDescriptors

import pandas as p
import os
path = os.path.join(os.getcwd(), 'C:\\Users\\Elif\\Desktop\\testReactionTable.csv')
df = p.read_csv(path)

def Convert(string): 
    list1=[] 
    list1[:0]=string 
    return list1 


def CreateProduct(stoichiometry):
    num=""
    targetFormula=""
    splitedFormula =""
    splitedFormula = stoichiometry.split("H")
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
                
    number= int(num)+2           
    targetFormula = "{0}{1}{2}{3}".format(splitedFormula[0],"H",str(number),secondPart)
    return targetFormula

for index, row in df.iterrows():
    try:
        splitedFormula = row['reactantStoichiometry'].split("H")
        if len(splitedFormula) ==1:
            createdReactantStoichiometry = row['productStoichiometry'].replace("H2", "")
            if row['reactantStoichiometry'] != createdReactantStoichiometry:
                print("*****************")
                print(index)
                print(row['reactantStoichiometry'])
                print(row['productStoichiometry'])
                print(createdProductStoichiometry)
                print("*****************")
        elif len(splitedFormula) ==2:
            createdProductStoichiometry = CreateProduct(row['reactantStoichiometry'])
            if row['productStoichiometry'] != createdProductStoichiometry:
                print("*****************")
                print(index)
                print(row['reactantStoichiometry'])
                print(row['productStoichiometry'])
                print(createdProductStoichiometry)
                print("*****************")
        else:
            print("BIG PROBLEM")
    except Exception as e:
        hey =1
        # print("--------------------------------")
        # print(index)
        # # print(row)
        # # print(e)
        # print("--------------------------------")
