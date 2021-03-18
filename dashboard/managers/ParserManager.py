import os
from .PathManager import ParserPaths
from .ObjectManager import FileFormat
from ..enums import RegexNames


import xml.etree.ElementTree as xml
import re as regx
import json

#https://www.freeformatter.com/xml-escape.html

class ParserManager():
                
    def Parse(commands, text):
        store=[]
        for regex in commands:
            tempDic=[]
            activeText=text
            for step in regex:
                try:
                    keys=step.find("keys")
                    command=regx.compile(step.find("command").text)
                    matchResult=command.finditer(activeText)
                    if regex[-1]==step:
                        for m in matchResult:
                            result=m.groupdict(keys.text)
                            store.append(result)
                    else:
                        for match in matchResult:
                            tempDic.append(match.groupdict(keys.text))
                        activeText=json.dumps(tempDic)
                except:
                    print('PARSER ERROR')
                    print("Regex Command :")
                    print(step.find("command").text)
                    print("Active Text :")
                    print(activeText)
                    continue        
        return store
        
    def GetRegexCommands(RegexName):
        path=ParserPaths.GetParserXml(RegexName)
        if os.path.exists(path):
            try:
                parsedXml=xml.parse(path)
                roots=parsedXml.getroot()
                return roots
            except:
                print('Could not parsed the regex file.')
        else:
            print('This path is not exist : ', path)
            
    def ResolveRegexName(calculationType):
        for key, value in RegexNames.__dict__.items():
            if key == calculationType:
                return value
            
