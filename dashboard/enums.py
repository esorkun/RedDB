# -*- coding: utf-8 -*-

#- regex ----------------------------------------------------------------------
class RegexNames:
    SinglePoint = 'regex-single'
    Optimization = 'regex-optimization'
    Smiles = 'regex-smiles'
    Pair = 'regex-pair'
    Configuration = 'config' 
    FileFormat = 'file-format'
    PairFileFormat = 'pair-file-format'

#- operating system -----------------------------------------------------------   
class FileExtensions:
    Output = '.out'
    Smiles = '.smi'
    Text = '.txt'
    
#- output file types ----------------------------------------------------------     
class FileTypes:
    SinglePoint = 'Single'
    Optimization= 'Optimization'
    Smiles = 'Smiles'
    
#- necessery file types -------------------------------------------------------
class AllowedFiles:
    smiles = 'smi'
    maestro= 'out'
    
#-  file types to read --------------------------------------------------------
class ToReadList:
    SinglePoint = 'Single'
    Optimization= 'Optimization'
    
