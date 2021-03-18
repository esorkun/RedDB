
import os    
        
class MoleculePaths:
       unzipedMoleculesPath = os.path.join(os.getcwd(), 'molecules')
       def ParentMoleculePath(parentName):
           return os.path.join(MoleculePaths.unzipedMoleculesPath, parentName)

class ParserPaths:
    configurationsPath = os.path.join(os.getcwd(), 'dashboard\\configurations')
    def GetParserXml(regexName):
        fileName='{}.xml'.format(regexName)
        return os.path.join(ParserPaths.configurationsPath, fileName)
    
class PairSmilesPaths:
       unzipedSmilesPath = os.path.join(os.getcwd(), 'pair')
       def ParentMoleculePath(parentName):
           return os.path.join(PairSmilesPaths.unzipedSmilesPath, parentName)
    
      
