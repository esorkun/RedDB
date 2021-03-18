# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 02:17:29 2020

@author: Elif
"""
import zipfile

class FileManager():
    def UnZipFiles(zipFilePath, extractTo):
      with zipfile.ZipFile(zipFilePath, 'r') as zip_ref:
          zip_ref.extractall(extractTo)

    def GetZippedFilesNames(zipFilePath):
        with zipfile.ZipFile(zipFilePath, 'r') as zip_ref:
            return zip_ref.namelist()

    def ReadFile(directoryPath):
      file = open(directoryPath, 'r')
      return file.read()

    def GetJobFileInfo(path):
      pass

