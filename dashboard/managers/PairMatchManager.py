# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 16:43:47 2020

@author: Elif
"""

try:
   from FileManager import FileManager
   from ObjectManager import MoleculeInfo, PairSmilesFileFormat
   from ParserManager import ParserManager
except ModuleNotFoundError:
   from .FileManager import FileManager
   from .ObjectManager import PairInfo, PairSmilesFileFormat   
   from .ParserManager import ParserManager

from .PathManager import PairSmilesPaths    
from ..models import PairMatchSmilesPackage
from ..enums import FileTypes, RegexNames, AllowedFiles, ToReadList
import os
   
class PairMatchManager():
    
    def PreparePackage(packageObj):
        resolvedFiles=[]
        fileNames = PairMatchManager.UnzipPackage(packageObj)
        parentMoleculeNames = PairMatchManager.GetParentNames(fileNames)        
        for parentName in parentMoleculeNames:
            moleculeFilePath = PairSmilesPaths.ParentMoleculePath(parentName)
            regexCommands = ParserManager.GetRegexCommands(RegexNames.PairFileFormat)
            resolvedFile = PairMatchManager.ResolveFileNames(parentName, moleculeFilePath, regexCommands)
            resolvedFiles.append(resolvedFile)
        return resolvedFiles

    def UnzipPackage(packageObj):
        dataPackage = PairMatchSmilesPackage.objects.get(id=packageObj.id)
        dataPackagePath = os.path.abspath(dataPackage.package.url)
        fileNames = FileManager.GetZippedFilesNames(dataPackagePath)
        # FileManager.UnZipFiles(dataPackagePath, PairSmilesPaths.unzipedSmilesPath)
        return fileNames
    
    def GetParentNames(filenames):
        mainMolecules =[]
        for file in filenames:
            parent = file.split('/')[0]
            mainMolecules.append(parent)
        return set(mainMolecules)
    
    def ResolveFileNames(parentName, moleculeFilePath, regexCommands): 
        pairSmileInfos = []
        for dirpath, foldernames, filenames in os.walk(moleculeFilePath):
            for file in filenames:
                fullPath = os.path.join(dirpath, file)
                if os.path.exists(fullPath):
                    pm = ParserManager.Parse(regexCommands, fullPath)
                    if len(pm) != 0:
                        pairFileDic ={}
                        for i in pm:
                            pairFileDic = {**pairFileDic, **i}
                        pairSmilesFile = PairSmilesFileFormat(**pairFileDic)
                        pairSmileInfos.append(pairSmilesFile)        
        return pairSmileInfos
    
    def ReadPairSmilesFile(regexCommands, pairSmilesFileObj):      
        pairParseResults = ParserManager.Parse(regexCommands, FileManager.ReadFile(pairSmilesFileObj.smilesFilePath))
        pairInfos =[]
        for result in pairParseResults:
            pairInfoDic = {}
            pairInfoDic = {**pairSmilesFileObj.__dict__, **result}
            pInfo = PairInfo(**pairInfoDic)
            pairInfos.append(pInfo)
        return pairInfos
    
    def ExtractData(packageObj):
        reacitonInfos = PairMatchManager.PreparePackage(packageObj)
        smilesRegexCommands = ParserManager.GetRegexCommands(RegexNames.Pair)
        pairSmilesInfos =[]
        for parent in reacitonInfos:
            for rInfo in parent:
                pairInfos = PairMatchManager.ReadPairSmilesFile(smilesRegexCommands, rInfo)
                pairSmilesInfos = pairSmilesInfos + pairInfos
        return pairSmilesInfos
    