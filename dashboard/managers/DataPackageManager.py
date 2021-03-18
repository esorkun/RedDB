import os

try:
   from FileManager import FileManager
   from ObjectManager import MoleculeInfo, FileFormat
   from ParserManager import ParserManager
except ModuleNotFoundError:
   from .FileManager import FileManager
   from .ObjectManager import MoleculeInfo, FileFormat   
   from .ParserManager import ParserManager


from .PathManager import MoleculePaths
from ..enums import FileTypes, RegexNames, AllowedFiles, ToReadList
from ..models import DataPackage


class DataPackageManager():
    
    def PreparePackage(packageObj):
        resolvedFiles=[]
        fileNames = DataPackageManager.UnzipPackage(packageObj)
        parentMoleculeNames = DataPackageManager.GetParentNames(fileNames)        
        for parentName in parentMoleculeNames:
            moleculeFilePath = MoleculePaths.ParentMoleculePath(parentName)
            regexCommands = ParserManager.GetRegexCommands(RegexNames.FileFormat)
            resolvedFile = DataPackageManager.ResolveFileNames(parentName, moleculeFilePath, regexCommands)
            resolvedFiles.append(resolvedFile)
        return resolvedFiles

    def UnzipPackage(packageObj):
        dataPackage = DataPackage.objects.get(id=packageObj.id)
        dataPackagePath = os.path.abspath(dataPackage.package.url)
        fileNames = FileManager.GetZippedFilesNames(dataPackagePath)
        # FileManager.UnZipFiles(dataPackagePath, MoleculePaths.unzipedMoleculesPath)
        return fileNames

    def GetParentNames(filenames):
        mainMolecules =[]
        for file in filenames:
            parent = file.split('/')[0]
            mainMolecules.append(parent)
        return set(mainMolecules)
    
    def ResolveFileNames(parentName, moleculeFilePath, regexCommands): 
        moleculeFiles = []
        for dirpath, foldernames, filenames in os.walk(moleculeFilePath):
          for file in filenames:
              fullPath = os.path.join(dirpath, file)
              if os.path.exists(fullPath):
                  pm = ParserManager.Parse(regexCommands, fullPath)
                  pm[0]['filePath'] = fullPath  
                  molFileInfo = FileFormat(**pm[0])
                  moleculeFiles.append(molFileInfo)
        return moleculeFiles
    
    def RemoveUnnecessaryFiles(datapackageFileInfo):
        allowedFilesInfo = []
        for parent in datapackageFileInfo:
            allowedFilesInfoParent =[]            
            for file in parent:
                for key, value in AllowedFiles.__dict__.items():
                    if key.startswith( '__' ) : continue
                    if (file.extension == value):
                        allowedFilesInfoParent.append(file)
                        break
            allowedFilesInfo.append(allowedFilesInfoParent)
        return allowedFilesInfo
      
    def GetSmilesFiles(parentFolder):
        return [moleculeFile for moleculeFile in parentFolder if (moleculeFile.calcType == FileTypes.Smiles)]
      
    def GetRemainingList(mainList, removeList):
        return [item for item in mainList if item not in removeList]
        
    def GetSmilesAndNames(smilesRegexCommands, smilesFiles):
        moleculesInfo = []
        for moleculeFile in smilesFiles:
            parsedSmilesAndNames = ParserManager.Parse(smilesRegexCommands, FileManager.ReadFile(moleculeFile.filePath))
            moleculeInfo = moleculeFile.__dict__
            for item in parsedSmilesAndNames:
                moleculeDict = {**moleculeInfo, **item}
                molecule = MoleculeInfo(**moleculeDict)
                moleculesInfo.append(molecule)
        return moleculesInfo
    
    def ReadCalculationFile(regexCommands, molecule, moleculeFileInfo):
            parsedCalculationValues = ParserManager.Parse(regexCommands,FileManager.ReadFile(moleculeFileInfo.filePath))    
            return parsedCalculationValues
    
    def OrganizeCalculationResults(results):
        calculationResults = {'Atoms': []}
        for item in results:
            if "atom" in item: 
                calculationResults['Atoms'].append(item)
            else: 
                calculationResults.update(item)
        return calculationResults
        
    def ExtractData(packageObj):
        datapackageFileInfo = DataPackageManager.PreparePackage(packageObj)
        filteredFilesInfo = DataPackageManager.RemoveUnnecessaryFiles(datapackageFileInfo)

        moleculesInfos = []              
        smilesRegexCommands = ParserManager.GetRegexCommands(RegexNames.Smiles)
        for parentFile in filteredFilesInfo:
            smilesFiles =DataPackageManager.GetSmilesFiles(parentFile)
            molecules = DataPackageManager.GetSmilesAndNames(smilesRegexCommands, smilesFiles)
            calculationFiles = DataPackageManager.GetRemainingList(parentFile, smilesFiles) # it removes smiles documents from the list.

            for file in calculationFiles:             
                for key, value in ToReadList.__dict__.items():
                    for mol in molecules:                                       
                        if (mol.parent == file.parent 
                            and mol.parentGivenId == file.parentGivenId
                            and mol.h == file.h
                            and mol.funcGroup == file.funcGroup
                            and mol.actualName in file.fileName
                            and file.calcType == value):
                            if hasattr(mol, value):
                                lst =[]
                                lst = getattr(mol, value)
                                lst.append(file)
                                delattr(mol, value) 
                                setattr(mol,value, lst)
                            else:
                                lst=[]
                                lst.append(file)                                
                                setattr(mol,value, lst)
                                                       
            for key, value in ToReadList.__dict__.items():
                if key.startswith( '__' ) : continue            
                regexCommands = ParserManager.GetRegexCommands(getattr(RegexNames, key))
                for mol in molecules:
                    moleculeFileInfos = getattr(mol, value, None)
                    if moleculeFileInfos == None:
                        organizedCalculationResults = None
                    else:
                        organizedCalculationResults =[]
                        for moleculeFileInfo in moleculeFileInfos:
                            calculationResults = DataPackageManager.ReadCalculationFile(regexCommands, mol, moleculeFileInfo)
                            if calculationResults == None : continue
                            organizedItem = DataPackageManager.OrganizeCalculationResults(calculationResults)
                            organizedItem.update(moleculeFileInfo.__dict__)
                            organizedCalculationResults.append(organizedItem)
                    attributeName = "{0}Results".format(value)
                    
                    setattr(mol,attributeName, organizedCalculationResults)
            moleculesInfo = {}           
            moleculesInfo["smilesFiles"] = smilesFiles
            moleculesInfo["molecules"] = molecules
            moleculesInfo["calculationFiles"] = calculationFiles 
            moleculesInfos.append(moleculesInfo)     
            
        return moleculesInfos
    
            
