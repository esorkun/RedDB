from django.contrib import admin
from .managers.DataPackageManager import DataPackageManager 
from .managers.MoleculeManager import MoleculeManager
from .managers.PairMatchManager import PairMatchManager


from dashboard.models import Reaction, PairMatchSmilesPackage, DataPackage, Molecule, JobSetting, OtherInfo, MoleculeInfo, FunctionalGroup, Job, ScfCalc, ChCalc, AtomicProperties, Optimization, CpolarCalc, PbfCalc, OptimizationGeometry

from operator import attrgetter
from itertools import groupby

import sys, os, traceback
from decimal import Decimal


class DataPackageAdmin(admin.ModelAdmin):
    
    def save_opt(molecule, moleculeModel, obj, optimizationResults):
        
        jobSettingOptimization, jsCreatedOptimization = JobSetting.objects.get_or_create(
            basisSet = optimizationResults["basisSet"],
            netMoleculerCharge = optimizationResults["netMolecularCharge"],
            multiplicity = optimizationResults["multiplicity"],
            scfCalculation = optimizationResults["scfCalculation"],
            dft = optimizationResults["dft"],
            maxScfIterations = optimizationResults["maxScfIterations"],
            solvent = None,
            internalDielectric = None,
            continuumDielectric = None,
            solventProbe = None, 
            pbfVersion = None)
                   
        functionalGroup, functionalGroupCreated = FunctionalGroup.objects.get_or_create(
            stoichiometry = optimizationResults["funcGroup"])
        
        moleculeInfoOptimizaiton, moleculeInfoOptCreated = MoleculeInfo.objects.get_or_create(
            moleculerWeight = optimizationResults["molecularWeight"],
            stoichiometry = optimizationResults["stoichiometry"])
        
        jobOptimization, jobOptimizationCreated = Job.objects.get_or_create(
            user = obj.user,
            jobType = "Optimization",
            jobId = optimizationResults["jobId"],
            calcNumber = optimizationResults["calcNumber"],
            jobSetting =jobSettingOptimization,
            moleculeInfo = moleculeInfoOptimizaiton,
            functionalGroup = functionalGroup,
            name = optimizationResults["fileName"],
            path = optimizationResults["filePath"],
            reactionStep = optimizationResults["h"],
            molecule = moleculeModel,
            defaults={'dataPackage': obj})
            
        optimization, optimizationCreated = Optimization.objects.get_or_create(
            job = jobOptimization,
            convergence = optimizationResults["Convergence"],
            optGasEnergy = optimizationResults["GasPhaseEnergy"],
            optGasIterationNumber = optimizationResults["GasPhaseScfIter"],
            optGasHomo = optimizationResults["GasPhaseHomoEnergy"],
            optGasLumo = optimizationResults["GasPhaseLumoEnergy"])
     
        atomicOptCalc = optimizationResults["Atoms"]
        for ac in atomicOptCalc:
            optGeo, optGeoCreated = OptimizationGeometry.objects.get_or_create(
                job = jobOptimization,
                atom = ac["atom"],
                atomId = ac["atomId"])            
            updateOperation = OptimizationGeometry.objects.filter(id = optGeo.id).update(**ac)

    def save_single(molecule, moleculeModel, obj,singlePointResults):
        functionalGroup, functionalGroupCreated = FunctionalGroup.objects.get_or_create(
            stoichiometry = singlePointResults["funcGroup"])

        jobSetting, jsCreated = JobSetting.objects.get_or_create(
            basisSet = singlePointResults["basisSet"],
            netMoleculerCharge = singlePointResults["netMolecularCharge"],
            multiplicity = singlePointResults["multiplicity"],
            solvent = singlePointResults["solvent"],
            scfCalculation = singlePointResults["scfCalculation"],
            dft = singlePointResults["dft"],
            solvationEnergy = singlePointResults["solvationEnergy"],
            hyperPolEqu = singlePointResults["hyperPolEqu"],
            maxScfIterations = singlePointResults["maxScfIterations"],
            internalDielectric = singlePointResults["internalDielectric"],
            continuumDielectric = singlePointResults["continuumDielectric"],
            solventProbe = singlePointResults["solventProbe"],
            pbfVersion = singlePointResults["pbfVersion"])

        moleculeInfo, moleculeInfoCreated = MoleculeInfo.objects.get_or_create(
            moleculerWeight = singlePointResults["molecularWeight"],
            stoichiometry = singlePointResults["stoichiometry"])

        job, jobCreated = Job.objects.get_or_create(
            user = obj.user,
            jobType = "SinglePoint",
            jobId = singlePointResults["jobId"],
            calcNumber = singlePointResults["calcNumber"],
            jobSetting =jobSetting,
            moleculeInfo = moleculeInfo,
            functionalGroup = functionalGroup,
            name = singlePointResults["fileName"],
            path = singlePointResults["filePath"],
            reactionStep = singlePointResults["h"],
            molecule = moleculeModel,
            defaults={'dataPackage': obj})

        atomicCalc = singlePointResults["Atoms"]
        
        for ac in atomicCalc:
            atomicP, apCreated = AtomicProperties.objects.get_or_create(
                job = job,
                atom = ac["atom"],
                atomId = ac["atomId"])            
            updateOperation = AtomicProperties.objects.filter(id = atomicP.id).update(**ac)
        
        try:
            OtherInfoSingle, OtherInfoSingleCreated = OtherInfo.objects.get_or_create(
                job = job,
                nuclearRepulsionEnergy = Decimal(singlePointResults["nuclearRepulsionEnergy"]),
                pointGroupUsed = singlePointResults["pointGroupUsed"],
                molecularPointGroup = singlePointResults["molecularPointGroup"])
            if not OtherInfoSingleCreated:
                OtherInfoSingle = None
                print("JobID : {0} --------- ERROR : OtherInfoSingle".format(job.id))
        except:
            OtherInfoSingle = None
            print("JobID : {0} --------- ERROR : OtherInfoSingle".format(job.id))
        
        try:
            pbf, pbfCreated = PbfCalc.objects.get_or_create( 
                job = job,
                molecularSurface = Decimal(singlePointResults["molecularSurface"]),
                reactionFieldEnergy = Decimal(singlePointResults["reactionFieldEnergy"]),
                solventAccessSurface = Decimal(singlePointResults["solventAccessSurface"]),
                cavityEnergy = Decimal(singlePointResults["cavityEnergy"]))
            if not pbfCreated:
                pbf = None
                print("JobID : {0} --------- ERROR : pbf".format(job.id))
        except:
            pbf = None
            print("JobID : {0} --------- ERROR : pbf".format(job.id))

        try:
            cPolar, cPolarCreated = CpolarCalc.objects.get_or_create(
                job = job,
                alpha = Decimal(singlePointResults["alpha"]),
                dalpha = Decimal(singlePointResults["dalpha"]),
                beta = singlePointResults["beta"])
            if not cPolarCreated:
                cPolar = None
                print("JobID : {0} --------- ERROR : cPolar".format(job.id))
        except:
            cPolar = None
            print("JobID : {0} --------- ERROR : cPolar".format(job.id))

        try:
            scfCalcResults, scfCalcCreated = ScfCalc.objects.get_or_create(
                job = job,
                gasEnergy = Decimal(singlePointResults["GasPhaseEnergy"]),
                gasIterationNumber = singlePointResults["GasPhaseScfIter"],
                gasHomo = Decimal(singlePointResults["GasPhaseHomoEnergy"]),
                gasLumo = Decimal(singlePointResults["GasPhaseLumoEnergy"]),
                solutionEnergy = Decimal(singlePointResults["SolutionPhaseEnergy"]),
                solutionIterationNumber = singlePointResults["scfIteration"],
                solutionHomo = Decimal(singlePointResults["solutionHomoEnergy"]),
                solutionLumo = Decimal(singlePointResults["solutionLumoEnergy"]))
            if not scfCalcCreated:
                scfCalcResults = None
                print("JobID : {0} --------- ERROR : scfCalc".format(job.id))
        except:
            scfCalcResults = None
            print("JobID : {0} --------- ERROR : scfCalc".format(job.id))

            
        
        try:
            chCalc, chCalcCreated = ChCalc.objects.get_or_create(
                job = job,
                mqmwDipoleMomentsX = Decimal(singlePointResults["mqmwDipoleMomentsX"]),
                mqmwDipoleMomentsY = Decimal(singlePointResults["mqmwDipoleMomentsY"]),
                mqmwDipoleMomentsZ = Decimal(singlePointResults["mqmwDipoleMomentsZ"]),
                mqmwDipoleMomentsTot = Decimal(singlePointResults["mqmwDipoleMomentsTot"]),
                mqmwQuadrupoleMomentsXX = Decimal(singlePointResults["mqmwQuadrupoleMomentsXX"]),
                mqmwQuadrupoleMomentsYY = Decimal(singlePointResults["mqmwQuadrupoleMomentsYY"]),
                mqmwQuadrupoleMomentsZZ = Decimal(singlePointResults["mqmwQuadrupoleMomentsZZ"]),
                mqmwQuadrupoleMomentsXY = Decimal(singlePointResults["mqmwQuadrupoleMomentsXY"]),
                mqmwQuadrupoleMomentsXZ = Decimal(singlePointResults["mqmwQuadrupoleMomentsXZ"]),
                mqmwQuadrupoleMomentsYZ = Decimal(singlePointResults["mqmwQuadrupoleMomentsYZ"]),
                mqmwTracelessQuadrupoleXXYY = Decimal(singlePointResults["mqmwTracelessQuadrupoleXXYY"]),
                mqmwTracelessQuadrupol2ZZXXYY = Decimal(singlePointResults["mqmwTracelessQuadrupol2ZZXXYY"]),
                mqmwTracelessQuadrupoleXY = Decimal(singlePointResults["mqmwTracelessQuadrupoleXY"]),
                mqmwTracelessQuadrupoleXZ = Decimal(singlePointResults["mqmwTracelessQuadrupoleXZ"]),
                mqmwTracelessQuadrupoleYZ = Decimal(singlePointResults["mqmwTracelessQuadrupoleYZ"]),
                mqmwOctapoleMomentsXXX = Decimal(singlePointResults["mqmwOctapoleMomentsXXX"]),
                mqmwOctapoleMomentsYYY = Decimal(singlePointResults["mqmwOctapoleMomentsYYY"]),
                mqmwOctapoleMomentsZZZ = Decimal(singlePointResults["mqmwOctapoleMomentsZZZ"]),
                mqmwOctapoleMomentsXYY = Decimal(singlePointResults["mqmwOctapoleMomentsXYY"]),
                mqmwOctapoleMomentsXXY = Decimal(singlePointResults["mqmwOctapoleMomentsXXY"]),
                mqmwOctapoleMomentsXXZ = Decimal(singlePointResults["mqmwOctapoleMomentsXXZ"]),
                mqmwOctapoleMomentsXZZ = Decimal(singlePointResults["mqmwOctapoleMomentsXZZ"]),
                mqmwOctapoleMomentsYZZ = Decimal(singlePointResults["mqmwOctapoleMomentsYZZ"]),
                mqmwOctapoleMomentsYYZ = Decimal(singlePointResults["mqmwOctapoleMomentsYYZ"]),
                mqmwOctapoleMomentsXYZ = Decimal(singlePointResults["mqmwOctapoleMomentsXYZ"]),
                mqmwTracelessOctapoleXXX = Decimal(singlePointResults["mqmwTracelessOctapoleXXX"]),
                mqmwTracelessOctapoleYYY = Decimal(singlePointResults["mqmwTracelessOctapoleYYY"]),
                mqmwTracelessOctapoleZZZ = Decimal(singlePointResults["mqmwTracelessOctapoleZZZ"]),
                mqmwTracelessOctapoleXYZ = Decimal(singlePointResults["mqmwTracelessOctapoleXYZ"]),
                mqmwTracelessOctapoleXYYXZZ = Decimal(singlePointResults["mqmwTracelessOctapoleXYYXZZ"]),
                mqmwTracelessOctapoleXXYYZZ = Decimal(singlePointResults["mqmwTracelessOctapoleXXYYZZ"]),
                mqmwTracelessOctapoleXXZYYZ = Decimal(singlePointResults["mqmwTracelessOctapoleXXZYYZ"]),
                mqmwHexadecapoleMomentsXXXX = Decimal(singlePointResults["mqmwHexadecapoleMomentsXXXX"]),
                mqmwHexadecapoleMomentsYYYY = Decimal(singlePointResults["mqmwHexadecapoleMomentsYYYY"]),
                mqmwHexadecapoleMomentsZZZZ = Decimal(singlePointResults["mqmwHexadecapoleMomentsZZZZ"]),
                mqmwHexadecapoleMomentsXXXY = Decimal(singlePointResults["mqmwHexadecapoleMomentsXXXY"]),
                mqmwHexadecapoleMomentsXXXZ = Decimal(singlePointResults["mqmwHexadecapoleMomentsXXXZ"]),
                mqmwHexadecapoleMomentsYYYX = Decimal(singlePointResults["mqmwHexadecapoleMomentsYYYX"]),
                mqmwHexadecapoleMomentsYYYZ = Decimal(singlePointResults["mqmwHexadecapoleMomentsYYYZ"]),
                mqmwHexadecapoleMomentsZZZX = Decimal(singlePointResults["mqmwHexadecapoleMomentsZZZX"]),
                mqmwHexadecapoleMomentsZZZY = Decimal(singlePointResults["mqmwHexadecapoleMomentsZZZY"]),
                mqmwHexadecapoleMomentsXXYY = Decimal(singlePointResults["mqmwHexadecapoleMomentsXXYY"]),
                mqmwHexadecapoleMomentsXXZZ = Decimal(singlePointResults["mqmwHexadecapoleMomentsXXZZ"]),
                mqmwHexadecapoleMomentsYYZZ = Decimal(singlePointResults["mqmwHexadecapoleMomentsYYZZ"]),
                mqmwHexadecapoleMomentsXXYZ = Decimal(singlePointResults["mqmwHexadecapoleMomentsXXYZ"]),
                mqmwHexadecapoleMomentsYYXZ = Decimal(singlePointResults["mqmwHexadecapoleMomentsYYXZ"]),
                mqmwHexadecapoleMomentsZZXY = Decimal(singlePointResults["mqmwHexadecapoleMomentsZZXY"]),
                mqmwGridpointsChargeFit = Decimal(singlePointResults["mqmwGridpointsChargeFit"]),
                mqmwPossibleMaximum = Decimal(singlePointResults["mqmwPossibleMaximum"]),
                
                mqmwGasDipoleMomentsX = Decimal(singlePointResults["mqmwGasDipoleMomentsX"]),
                mqmwGasDipoleMomentsY = Decimal(singlePointResults["mqmwGasDipoleMomentsY"]),
                mqmwGasDipoleMomentsZ = Decimal(singlePointResults["mqmwGasDipoleMomentsZ"]),
                mqmwGasDipoleMomentsTot = Decimal(singlePointResults["mqmwGasDipoleMomentsTot"]),
                mqmwGasQuadrupoleMomentsXX = Decimal(singlePointResults["mqmwGasQuadrupoleMomentsXX"]),
                mqmwGasQuadrupoleMomentsYY = Decimal(singlePointResults["mqmwGasQuadrupoleMomentsYY"]),
                mqmwGasQuadrupoleMomentsZZ = Decimal(singlePointResults["mqmwGasQuadrupoleMomentsZZ"]),
                mqmwGasQuadrupoleMomentsXY = Decimal(singlePointResults["mqmwGasQuadrupoleMomentsXY"]),
                mqmwGasQuadrupoleMomentsXZ = Decimal(singlePointResults["mqmwGasQuadrupoleMomentsXZ"]),
                mqmwGasQuadrupoleMomentsYZ = Decimal(singlePointResults["mqmwGasQuadrupoleMomentsYZ"]),
                mqmwGasTracelessQuadrupoleXXYY = Decimal(singlePointResults["mqmwGasTracelessQuadrupoleXXYY"]),
                mqmwGasTracelessQuadrupol2ZZXXYY = Decimal(singlePointResults["mqmwGasTracelessQuadrupol2ZZXXYY"]),
                mqmwGasTracelessQuadrupoleXY = Decimal(singlePointResults["mqmwGasTracelessQuadrupoleXY"]),
                mqmwGasTracelessQuadrupoleXZ = Decimal(singlePointResults["mqmwGasTracelessQuadrupoleXZ"]),
                mqmwGasTracelessQuadrupoleYZ = Decimal(singlePointResults["mqmwGasTracelessQuadrupoleYZ"]),
                mqmwGasOctapoleMomentsXXX = Decimal(singlePointResults["mqmwGasOctapoleMomentsXXX"]),
                mqmwGasOctapoleMomentsYYY = Decimal(singlePointResults["mqmwGasOctapoleMomentsYYY"]),
                mqmwGasOctapoleMomentsZZZ = Decimal(singlePointResults["mqmwGasOctapoleMomentsZZZ"]),
                mqmwGasOctapoleMomentsXYY = Decimal(singlePointResults["mqmwGasOctapoleMomentsXYY"]),
                mqmwGasOctapoleMomentsXXY = Decimal(singlePointResults["mqmwGasOctapoleMomentsXXY"]),
                mqmwGasOctapoleMomentsXXZ = Decimal(singlePointResults["mqmwGasOctapoleMomentsXXZ"]),
                mqmwGasOctapoleMomentsXZZ = Decimal(singlePointResults["mqmwGasOctapoleMomentsXZZ"]),
                mqmwGasOctapoleMomentsYZZ = Decimal(singlePointResults["mqmwGasOctapoleMomentsYZZ"]),
                mqmwGasOctapoleMomentsYYZ = Decimal(singlePointResults["mqmwGasOctapoleMomentsYYZ"]),
                mqmwGasOctapoleMomentsXYZ = Decimal(singlePointResults["mqmwGasOctapoleMomentsXYZ"]),
                mqmwGasTracelessOctapoleXXX = Decimal(singlePointResults["mqmwGasTracelessOctapoleXXX"]),
                mqmwGasTracelessOctapoleYYY = Decimal(singlePointResults["mqmwGasTracelessOctapoleYYY"]),
                mqmwGasTracelessOctapoleZZZ = Decimal(singlePointResults["mqmwGasTracelessOctapoleZZZ"]),
                mqmwGasTracelessOctapoleXYZ = Decimal(singlePointResults["mqmwGasTracelessOctapoleXYZ"]),
                mqmwGasTracelessOctapoleXYYXZZ = Decimal(singlePointResults["mqmwGasTracelessOctapoleXYYXZZ"]),
                mqmwGasTracelessOctapoleXXYYZZ = Decimal(singlePointResults["mqmwGasTracelessOctapoleXXYYZZ"]),
                mqmwGasTracelessOctapoleXXZYYZ = Decimal(singlePointResults["mqmwGasTracelessOctapoleXXZYYZ"]),
                mqmwGasHexadecapoleMomentsXXXX = Decimal(singlePointResults["mqmwGasHexadecapoleMomentsXXXX"]),
                mqmwGasHexadecapoleMomentsYYYY = Decimal(singlePointResults["mqmwGasHexadecapoleMomentsYYYY"]),
                mqmwGasHexadecapoleMomentsZZZZ = Decimal(singlePointResults["mqmwGasHexadecapoleMomentsZZZZ"]),
                mqmwGasHexadecapoleMomentsXXXY = Decimal(singlePointResults["mqmwGasHexadecapoleMomentsXXXY"]),
                mqmwGasHexadecapoleMomentsXXXZ = Decimal(singlePointResults["mqmwGasHexadecapoleMomentsXXXZ"]),
                mqmwGasHexadecapoleMomentsYYYX = Decimal(singlePointResults["mqmwGasHexadecapoleMomentsYYYX"]),
                mqmwGasHexadecapoleMomentsYYYZ = Decimal(singlePointResults["mqmwGasHexadecapoleMomentsYYYZ"]),
                mqmwGasHexadecapoleMomentsZZZX = Decimal(singlePointResults["mqmwGasHexadecapoleMomentsZZZX"]),
                mqmwGasHexadecapoleMomentsZZZY = Decimal(singlePointResults["mqmwGasHexadecapoleMomentsZZZY"]),
                mqmwGasHexadecapoleMomentsXXYY = Decimal(singlePointResults["mqmwGasHexadecapoleMomentsXXYY"]),
                mqmwGasHexadecapoleMomentsXXZZ = Decimal(singlePointResults["mqmwGasHexadecapoleMomentsXXZZ"]),
                mqmwGasHexadecapoleMomentsYYZZ = Decimal(singlePointResults["mqmwGasHexadecapoleMomentsYYZZ"]),
                mqmwGasHexadecapoleMomentsXXYZ = Decimal(singlePointResults["mqmwGasHexadecapoleMomentsXXYZ"]),
                mqmwGasHexadecapoleMomentsYYXZ = Decimal(singlePointResults["mqmwGasHexadecapoleMomentsYYXZ"]),
                mqmwGasHexadecapoleMomentsZZXY = Decimal(singlePointResults["mqmwGasHexadecapoleMomentsZZXY"]),
                mqmwGasGridpointsChargeFit = Decimal(singlePointResults["mqmwGasGridpointsChargeFit"]),
                mqmwGasPossibleMaximum = Decimal(singlePointResults["mqmwGasPossibleMaximum"]),
                
                mepcGasDipoleMomentsX = Decimal(singlePointResults["mepcGasDipoleMomentsX"]),
                mepcGasDipoleMomentsY = Decimal(singlePointResults["mepcGasDipoleMomentsY"]),
                mepcGasDipoleMomentsZ = Decimal(singlePointResults["mepcGasDipoleMomentsZ"]),
                mepcGasDipoleMomentsTot = Decimal(singlePointResults["mepcGasDipoleMomentsTot"]),
                mepcGasQuadrupoleMomentsXX = Decimal(singlePointResults["mepcGasQuadrupoleMomentsXX"]),
                mepcGasQuadrupoleMomentsYY = Decimal(singlePointResults["mepcGasQuadrupoleMomentsYY"]),
                mepcGasQuadrupoleMomentsZZ = Decimal(singlePointResults["mepcGasQuadrupoleMomentsZZ"]),
                mepcGasQuadrupoleMomentsXY = Decimal(singlePointResults["mepcGasQuadrupoleMomentsXY"]),
                mepcGasQuadrupoleMomentsXZ = Decimal(singlePointResults["mepcGasQuadrupoleMomentsXZ"]),
                mepcGasQuadrupoleMomentsYZ = Decimal(singlePointResults["mepcGasQuadrupoleMomentsYZ"]),
                mepcGasTracelessQuadrupoleXXYY = Decimal(singlePointResults["mepcGasTracelessQuadrupoleXXYY"]),
                mepcGasTracelessQuadrupol2ZZXXYY = Decimal(singlePointResults["mepcGasTracelessQuadrupol2ZZXXYY"]),
                mepcGasTracelessQuadrupoleXY = Decimal(singlePointResults["mepcGasTracelessQuadrupoleXY"]),
                mepcGasTracelessQuadrupoleXZ = Decimal(singlePointResults["mepcGasTracelessQuadrupoleXZ"]),
                mepcGasTracelessQuadrupoleYZ = Decimal(singlePointResults["mepcGasTracelessQuadrupoleYZ"]),
                mepcGasOctapoleMomentsXXX = Decimal(singlePointResults["mepcGasOctapoleMomentsXXX"]),
                mepcGasOctapoleMomentsYYY = Decimal(singlePointResults["mepcGasOctapoleMomentsYYY"]),
                mepcGasOctapoleMomentsZZZ = Decimal(singlePointResults["mepcGasOctapoleMomentsZZZ"]),
                mepcGasOctapoleMomentsXYY = Decimal(singlePointResults["mepcGasOctapoleMomentsXYY"]),
                mepcGasOctapoleMomentsXXY = Decimal(singlePointResults["mepcGasOctapoleMomentsXXY"]),
                mepcGasOctapoleMomentsXXZ = Decimal(singlePointResults["mepcGasOctapoleMomentsXXZ"]),
                mepcGasOctapoleMomentsXZZ = Decimal(singlePointResults["mepcGasOctapoleMomentsXZZ"]),
                mepcGasOctapoleMomentsYZZ = Decimal(singlePointResults["mepcGasOctapoleMomentsYZZ"]),
                mepcGasOctapoleMomentsYYZ = Decimal(singlePointResults["mepcGasOctapoleMomentsYYZ"]),
                mepcGasOctapoleMomentsXYZ = Decimal(singlePointResults["mepcGasOctapoleMomentsXYZ"]),
                mepcGasTracelessOctapoleXXX = Decimal(singlePointResults["mepcGasTracelessOctapoleXXX"]),
                mepcGasTracelessOctapoleYYY = Decimal(singlePointResults["mepcGasTracelessOctapoleYYY"]),
                mepcGasTracelessOctapoleZZZ = Decimal(singlePointResults["mepcGasTracelessOctapoleZZZ"]),
                mepcGasTracelessOctapoleXYZ = Decimal(singlePointResults["mepcGasTracelessOctapoleXYZ"]),
                mepcGasTracelessOctapoleXYYXZZ = Decimal(singlePointResults["mepcGasTracelessOctapoleXYYXZZ"]),
                mepcGasTracelessOctapoleXXYYZZ = Decimal(singlePointResults["mepcGasTracelessOctapoleXXYYZZ"]),
                mepcGasTracelessOctapoleXXZYYZ = Decimal(singlePointResults["mepcGasTracelessOctapoleXXZYYZ"]),
                mepcGasHexadecapoleMomentsXXXX = Decimal(singlePointResults["mepcGasHexadecapoleMomentsXXXX"]),
                mepcGasHexadecapoleMomentsYYYY = Decimal(singlePointResults["mepcGasHexadecapoleMomentsYYYY"]),
                mepcGasHexadecapoleMomentsZZZZ = Decimal(singlePointResults["mepcGasHexadecapoleMomentsZZZZ"]),
                mepcGasHexadecapoleMomentsXXXY = Decimal(singlePointResults["mepcGasHexadecapoleMomentsXXXY"]),
                mepcGasHexadecapoleMomentsXXXZ = Decimal(singlePointResults["mepcGasHexadecapoleMomentsXXXZ"]),
                mepcGasHexadecapoleMomentsYYYX = Decimal(singlePointResults["mepcGasHexadecapoleMomentsYYYX"]),
                mepcGasHexadecapoleMomentsYYYZ = Decimal(singlePointResults["mepcGasHexadecapoleMomentsYYYZ"]),
                mepcGasHexadecapoleMomentsZZZX = Decimal(singlePointResults["mepcGasHexadecapoleMomentsZZZX"]),
                mepcGasHexadecapoleMomentsZZZY = Decimal(singlePointResults["mepcGasHexadecapoleMomentsZZZY"]),
                mepcGasHexadecapoleMomentsXXYY = Decimal(singlePointResults["mepcGasHexadecapoleMomentsXXYY"]),
                mepcGasHexadecapoleMomentsXXZZ = Decimal(singlePointResults["mepcGasHexadecapoleMomentsXXZZ"]),
                mepcGasHexadecapoleMomentsYYZZ = Decimal(singlePointResults["mepcGasHexadecapoleMomentsYYZZ"]),
                mepcGasHexadecapoleMomentsXXYZ = Decimal(singlePointResults["mepcGasHexadecapoleMomentsXXYZ"]),
                mepcGasHexadecapoleMomentsYYXZ = Decimal(singlePointResults["mepcGasHexadecapoleMomentsYYXZ"]),
                mepcGasHexadecapoleMomentsZZXY = Decimal(singlePointResults["mepcGasHexadecapoleMomentsZZXY"]),
                
                mepcDipoleMomentsX = Decimal(singlePointResults["mepcDipoleMomentsX"]),
                mepcDipoleMomentsY = Decimal(singlePointResults["mepcDipoleMomentsY"]),
                mepcDipoleMomentsZ = Decimal(singlePointResults["mepcDipoleMomentsZ"]),
                mepcDipoleMomentsTot = Decimal(singlePointResults["mepcDipoleMomentsTot"]),
                mepcQuadrupoleMomentsXX = Decimal(singlePointResults["mepcQuadrupoleMomentsXX"]),
                mepcQuadrupoleMomentsYY = Decimal(singlePointResults["mepcQuadrupoleMomentsYY"]),
                mepcQuadrupoleMomentsZZ = Decimal(singlePointResults["mepcQuadrupoleMomentsZZ"]),
                mepcQuadrupoleMomentsXY = Decimal(singlePointResults["mepcQuadrupoleMomentsXY"]),
                mepcQuadrupoleMomentsXZ = Decimal(singlePointResults["mepcQuadrupoleMomentsXZ"]),
                mepcQuadrupoleMomentsYZ = Decimal(singlePointResults["mepcQuadrupoleMomentsYZ"]),
                mepcTracelessQuadrupoleXXYY = Decimal(singlePointResults["mepcTracelessQuadrupoleXXYY"]),
                mepcTracelessQuadrupol2ZZXXYY = Decimal(singlePointResults["mepcTracelessQuadrupol2ZZXXYY"]),
                mepcTracelessQuadrupoleXY = Decimal(singlePointResults["mepcTracelessQuadrupoleXY"]),
                mepcTracelessQuadrupoleXZ = Decimal(singlePointResults["mepcTracelessQuadrupoleXZ"]),
                mepcTracelessQuadrupoleYZ = Decimal(singlePointResults["mepcTracelessQuadrupoleYZ"]),
                mepcOctapoleMomentsXXX = Decimal(singlePointResults["mepcOctapoleMomentsXXX"]),
                mepcOctapoleMomentsYYY = Decimal(singlePointResults["mepcOctapoleMomentsYYY"]),
                mepcOctapoleMomentsZZZ = Decimal(singlePointResults["mepcOctapoleMomentsZZZ"]),
                mepcOctapoleMomentsXYY = Decimal(singlePointResults["mepcOctapoleMomentsXYY"]),
                mepcOctapoleMomentsXXY = Decimal(singlePointResults["mepcOctapoleMomentsXXY"]),
                mepcOctapoleMomentsXXZ = Decimal(singlePointResults["mepcOctapoleMomentsXXZ"]),
                mepcOctapoleMomentsXZZ = Decimal(singlePointResults["mepcOctapoleMomentsXZZ"]),
                mepcOctapoleMomentsYZZ = Decimal(singlePointResults["mepcOctapoleMomentsYZZ"]),
                mepcOctapoleMomentsYYZ = Decimal(singlePointResults["mepcOctapoleMomentsYYZ"]),
                mepcOctapoleMomentsXYZ = Decimal(singlePointResults["mepcOctapoleMomentsXYZ"]),
                mepcTracelessOctapoleXXX = Decimal(singlePointResults["mepcTracelessOctapoleXXX"]),
                mepcTracelessOctapoleYYY = Decimal(singlePointResults["mepcTracelessOctapoleYYY"]),
                mepcTracelessOctapoleZZZ = Decimal(singlePointResults["mepcTracelessOctapoleZZZ"]),
                mepcTracelessOctapoleXYZ = Decimal(singlePointResults["mepcTracelessOctapoleXYZ"]),
                mepcTracelessOctapoleXYYXZZ = Decimal(singlePointResults["mepcTracelessOctapoleXYYXZZ"]),
                mepcTracelessOctapoleXXYYZZ = Decimal(singlePointResults["mepcTracelessOctapoleXXYYZZ"]),
                mepcTracelessOctapoleXXZYYZ = Decimal(singlePointResults["mepcTracelessOctapoleXXZYYZ"]),
                mepcHexadecapoleMomentsXXXX = Decimal(singlePointResults["mepcHexadecapoleMomentsXXXX"]),
                mepcHexadecapoleMomentsYYYY = Decimal(singlePointResults["mepcHexadecapoleMomentsYYYY"]),
                mepcHexadecapoleMomentsZZZZ = Decimal(singlePointResults["mepcHexadecapoleMomentsZZZZ"]),
                mepcHexadecapoleMomentsXXXY = Decimal(singlePointResults["mepcHexadecapoleMomentsXXXY"]),
                mepcHexadecapoleMomentsXXXZ = Decimal(singlePointResults["mepcHexadecapoleMomentsXXXZ"]),
                mepcHexadecapoleMomentsYYYX = Decimal(singlePointResults["mepcHexadecapoleMomentsYYYX"]),
                mepcHexadecapoleMomentsYYYZ = Decimal(singlePointResults["mepcHexadecapoleMomentsYYYZ"]),
                mepcHexadecapoleMomentsZZZX = Decimal(singlePointResults["mepcHexadecapoleMomentsZZZX"]),
                mepcHexadecapoleMomentsZZZY = Decimal(singlePointResults["mepcHexadecapoleMomentsZZZY"]),
                mepcHexadecapoleMomentsXXYY = Decimal(singlePointResults["mepcHexadecapoleMomentsXXYY"]),
                mepcHexadecapoleMomentsXXZZ = Decimal(singlePointResults["mepcHexadecapoleMomentsXXZZ"]),
                mepcHexadecapoleMomentsYYZZ = Decimal(singlePointResults["mepcHexadecapoleMomentsYYZZ"]),
                mepcHexadecapoleMomentsXXYZ = Decimal(singlePointResults["mepcHexadecapoleMomentsXXYZ"]),
                mepcHexadecapoleMomentsYYXZ = Decimal(singlePointResults["mepcHexadecapoleMomentsYYXZ"]),
                mepcHexadecapoleMomentsZZXY = Decimal(singlePointResults["mepcHexadecapoleMomentsZZXY"]),
                
                mmcGasDipoleMomentsX = Decimal(singlePointResults["mmcGasDipoleMomentsX"]),
                mmcGasDipoleMomentsY = Decimal(singlePointResults["mmcGasDipoleMomentsY"]),
                mmcGasDipoleMomentsZ = Decimal(singlePointResults["mmcGasDipoleMomentsZ"]),
                mmcGasDipoleMomentsTot = Decimal(singlePointResults["mmcGasDipoleMomentsTot"]),
                mmcGasQuadrupoleMomentsXX = Decimal(singlePointResults["mmcGasQuadrupoleMomentsXX"]),
                mmcGasQuadrupoleMomentsYY = Decimal(singlePointResults["mmcGasQuadrupoleMomentsYY"]),
                mmcGasQuadrupoleMomentsZZ = Decimal(singlePointResults["mmcGasQuadrupoleMomentsZZ"]),
                mmcGasQuadrupoleMomentsXY = Decimal(singlePointResults["mmcGasQuadrupoleMomentsXY"]),
                mmcGasQuadrupoleMomentsXZ = Decimal(singlePointResults["mmcGasQuadrupoleMomentsXZ"]),
                mmcGasQuadrupoleMomentsYZ = Decimal(singlePointResults["mmcGasQuadrupoleMomentsYZ"]),
                mmcGasTracelessQuadrupoleXXYY = Decimal(singlePointResults["mmcGasTracelessQuadrupoleXXYY"]),
                mmcGasTracelessQuadrupol2ZZXXYY = Decimal(singlePointResults["mmcGasTracelessQuadrupol2ZZXXYY"]),
                mmcGasTracelessQuadrupoleXY = Decimal(singlePointResults["mmcGasTracelessQuadrupoleXY"]),
                mmcGasTracelessQuadrupoleXZ = Decimal(singlePointResults["mmcGasTracelessQuadrupoleXZ"]),
                mmcGasTracelessQuadrupoleYZ = Decimal(singlePointResults["mmcGasTracelessQuadrupoleYZ"]),
                mmcGasOctapoleMomentsXXX = Decimal(singlePointResults["mmcGasOctapoleMomentsXXX"]),
                mmcGasOctapoleMomentsYYY = Decimal(singlePointResults["mmcGasOctapoleMomentsYYY"]),
                mmcGasOctapoleMomentsZZZ = Decimal(singlePointResults["mmcGasOctapoleMomentsZZZ"]),
                mmcGasOctapoleMomentsXYY = Decimal(singlePointResults["mmcGasOctapoleMomentsXYY"]),
                mmcGasOctapoleMomentsXXY = Decimal(singlePointResults["mmcGasOctapoleMomentsXXY"]),
                mmcGasOctapoleMomentsXXZ = Decimal(singlePointResults["mmcGasOctapoleMomentsXXZ"]),
                mmcGasOctapoleMomentsXZZ = Decimal(singlePointResults["mmcGasOctapoleMomentsXZZ"]),
                mmcGasOctapoleMomentsYZZ = Decimal(singlePointResults["mmcGasOctapoleMomentsYZZ"]),
                mmcGasOctapoleMomentsYYZ = Decimal(singlePointResults["mmcGasOctapoleMomentsYYZ"]),
                mmcGasOctapoleMomentsXYZ = Decimal(singlePointResults["mmcGasOctapoleMomentsXYZ"]),
                mmcGasTracelessOctapoleXXX = Decimal(singlePointResults["mmcGasTracelessOctapoleXXX"]),
                mmcGasTracelessOctapoleYYY = Decimal(singlePointResults["mmcGasTracelessOctapoleYYY"]),
                mmcGasTracelessOctapoleZZZ = Decimal(singlePointResults["mmcGasTracelessOctapoleZZZ"]),
                mmcGasTracelessOctapoleXYZ = Decimal(singlePointResults["mmcGasTracelessOctapoleXYZ"]),
                mmcGasTracelessOctapoleXYYXZZ = Decimal(singlePointResults["mmcGasTracelessOctapoleXYYXZZ"]),
                mmcGasTracelessOctapoleXXYYZZ = Decimal(singlePointResults["mmcGasTracelessOctapoleXXYYZZ"]),
                mmcGasTracelessOctapoleXXZYYZ = Decimal(singlePointResults["mmcGasTracelessOctapoleXXZYYZ"]),
                mmcGasHexadecapoleMomentsXXXX = Decimal(singlePointResults["mmcGasHexadecapoleMomentsXXXX"]),
                mmcGasHexadecapoleMomentsYYYY = Decimal(singlePointResults["mmcGasHexadecapoleMomentsYYYY"]),
                mmcGasHexadecapoleMomentsZZZZ = Decimal(singlePointResults["mmcGasHexadecapoleMomentsZZZZ"]),
                mmcGasHexadecapoleMomentsXXXY = Decimal(singlePointResults["mmcGasHexadecapoleMomentsXXXY"]),
                mmcGasHexadecapoleMomentsXXXZ = Decimal(singlePointResults["mmcGasHexadecapoleMomentsXXXZ"]),
                mmcGasHexadecapoleMomentsYYYX = Decimal(singlePointResults["mmcGasHexadecapoleMomentsYYYX"]),
                mmcGasHexadecapoleMomentsYYYZ = Decimal(singlePointResults["mmcGasHexadecapoleMomentsYYYZ"]),
                mmcGasHexadecapoleMomentsZZZX = Decimal(singlePointResults["mmcGasHexadecapoleMomentsZZZX"]),
                mmcGasHexadecapoleMomentsZZZY = Decimal(singlePointResults["mmcGasHexadecapoleMomentsZZZY"]),
                mmcGasHexadecapoleMomentsXXYY = Decimal(singlePointResults["mmcGasHexadecapoleMomentsXXYY"]),
                mmcGasHexadecapoleMomentsXXZZ = Decimal(singlePointResults["mmcGasHexadecapoleMomentsXXZZ"]),
                mmcGasHexadecapoleMomentsYYZZ = Decimal(singlePointResults["mmcGasHexadecapoleMomentsYYZZ"]),
                mmcGasHexadecapoleMomentsXXYZ = Decimal(singlePointResults["mmcGasHexadecapoleMomentsXXYZ"]),
                mmcGasHexadecapoleMomentsYYXZ = Decimal(singlePointResults["mmcGasHexadecapoleMomentsYYXZ"]),
                mmcGasHexadecapoleMomentsZZXY = Decimal(singlePointResults["mmcGasHexadecapoleMomentsZZXY"]),
                
                mmcDipoleMomentsX = Decimal(singlePointResults["mmcDipoleMomentsX"]),
                mmcDipoleMomentsY = Decimal(singlePointResults["mmcDipoleMomentsY"]),
                mmcDipoleMomentsZ = Decimal(singlePointResults["mmcDipoleMomentsZ"]),
                mmcDipoleMomentsTot = Decimal(singlePointResults["mmcDipoleMomentsTot"]),
                mmcQuadrupoleMomentsXX = Decimal(singlePointResults["mmcQuadrupoleMomentsXX"]),
                mmcQuadrupoleMomentsYY = Decimal(singlePointResults["mmcQuadrupoleMomentsYY"]),
                mmcQuadrupoleMomentsZZ = Decimal(singlePointResults["mmcQuadrupoleMomentsZZ"]),
                mmcQuadrupoleMomentsXY = Decimal(singlePointResults["mmcQuadrupoleMomentsXY"]),
                mmcQuadrupoleMomentsXZ = Decimal(singlePointResults["mmcQuadrupoleMomentsXZ"]),
                mmcQuadrupoleMomentsYZ = Decimal(singlePointResults["mmcQuadrupoleMomentsYZ"]),
                mmcTracelessQuadrupoleXXYY = Decimal(singlePointResults["mmcTracelessQuadrupoleXXYY"]),
                mmcTracelessQuadrupol2ZZXXYY = Decimal(singlePointResults["mmcTracelessQuadrupol2ZZXXYY"]),
                mmcTracelessQuadrupoleXY = Decimal(singlePointResults["mmcTracelessQuadrupoleXY"]),
                mmcTracelessQuadrupoleXZ = Decimal(singlePointResults["mmcTracelessQuadrupoleXZ"]),
                mmcTracelessQuadrupoleYZ = Decimal(singlePointResults["mmcTracelessQuadrupoleYZ"]),
                mmcOctapoleMomentsXXX = Decimal(singlePointResults["mmcOctapoleMomentsXXX"]),
                mmcOctapoleMomentsYYY = Decimal(singlePointResults["mmcOctapoleMomentsYYY"]),
                mmcOctapoleMomentsZZZ = Decimal(singlePointResults["mmcOctapoleMomentsZZZ"]),
                mmcOctapoleMomentsXYY = Decimal(singlePointResults["mmcOctapoleMomentsXYY"]),
                mmcOctapoleMomentsXXY = Decimal(singlePointResults["mmcOctapoleMomentsXXY"]),
                mmcOctapoleMomentsXXZ = Decimal(singlePointResults["mmcOctapoleMomentsXXZ"]),
                mmcOctapoleMomentsXZZ = Decimal(singlePointResults["mmcOctapoleMomentsXZZ"]),
                mmcOctapoleMomentsYZZ = Decimal(singlePointResults["mmcOctapoleMomentsYZZ"]),
                mmcOctapoleMomentsYYZ = Decimal(singlePointResults["mmcOctapoleMomentsYYZ"]),
                mmcOctapoleMomentsXYZ = Decimal(singlePointResults["mmcOctapoleMomentsXYZ"]),
                mmcTracelessOctapoleXXX = Decimal(singlePointResults["mmcTracelessOctapoleXXX"]),
                mmcTracelessOctapoleYYY = Decimal(singlePointResults["mmcTracelessOctapoleYYY"]),
                mmcTracelessOctapoleZZZ = Decimal(singlePointResults["mmcTracelessOctapoleZZZ"]),
                mmcTracelessOctapoleXYZ = Decimal(singlePointResults["mmcTracelessOctapoleXYZ"]),
                mmcTracelessOctapoleXYYXZZ = Decimal(singlePointResults["mmcTracelessOctapoleXYYXZZ"]),
                mmcTracelessOctapoleXXYYZZ = Decimal(singlePointResults["mmcTracelessOctapoleXXYYZZ"]),
                mmcTracelessOctapoleXXZYYZ = Decimal(singlePointResults["mmcTracelessOctapoleXXZYYZ"]),
                mmcHexadecapoleMomentsXXXX = Decimal(singlePointResults["mmcHexadecapoleMomentsXXXX"]),
                mmcHexadecapoleMomentsYYYY = Decimal(singlePointResults["mmcHexadecapoleMomentsYYYY"]),
                mmcHexadecapoleMomentsZZZZ = Decimal(singlePointResults["mmcHexadecapoleMomentsZZZZ"]),
                mmcHexadecapoleMomentsXXXY = Decimal(singlePointResults["mmcHexadecapoleMomentsXXXY"]),
                mmcHexadecapoleMomentsXXXZ = Decimal(singlePointResults["mmcHexadecapoleMomentsXXXZ"]),
                mmcHexadecapoleMomentsYYYX = Decimal(singlePointResults["mmcHexadecapoleMomentsYYYX"]),
                mmcHexadecapoleMomentsYYYZ = Decimal(singlePointResults["mmcHexadecapoleMomentsYYYZ"]),
                mmcHexadecapoleMomentsZZZX = Decimal(singlePointResults["mmcHexadecapoleMomentsZZZX"]),
                mmcHexadecapoleMomentsZZZY = Decimal(singlePointResults["mmcHexadecapoleMomentsZZZY"]),
                mmcHexadecapoleMomentsXXYY = Decimal(singlePointResults["mmcHexadecapoleMomentsXXYY"]),
                mmcHexadecapoleMomentsXXZZ = Decimal(singlePointResults["mmcHexadecapoleMomentsXXZZ"]),
                mmcHexadecapoleMomentsYYZZ = Decimal(singlePointResults["mmcHexadecapoleMomentsYYZZ"]),
                mmcHexadecapoleMomentsXXYZ = Decimal(singlePointResults["mmcHexadecapoleMomentsXXYZ"]),
                mmcHexadecapoleMomentsYYXZ = Decimal(singlePointResults["mmcHexadecapoleMomentsYYXZ"]),
                mmcHexadecapoleMomentsZZXY = Decimal(singlePointResults["mmcHexadecapoleMomentsZZXY"]))
            if not chCalcCreated:
                chCalc = None
                print("JobID : {0} --------- ERROR : chCalc".format(job.id))
        except:
            chCalc = None
            print("JobID : {0} --------- ERROR : chCalc".format(job.id))
        
            
        
    def save_calc(molecule, moleculeModel, obj):
        if hasattr(molecule, 'OptimizationResults'):
            optimizationResults = getattr(molecule, 'OptimizationResults')
            if optimizationResults != None:
                for item in optimizationResults:
                    DataPackageAdmin.save_opt(molecule, moleculeModel, obj, item)
        else:
            print("Missing OptimizationResults : ")
            print(molecule.__dict__)
         
        if hasattr(molecule, 'SingleResults'):
            singlePointResults = getattr(molecule, 'SingleResults')
            if singlePointResults != None:
                for item in singlePointResults:
                    DataPackageAdmin.save_single(molecule, moleculeModel, obj, item)
        else:
            print("Missing SingleResults : ")
            print(molecule.__dict__)
       

    def save_model(self, request, obj, form, change):
        obj.user = request.user
        if (True): #ValidationFilesManager.ValidateData(obj)
            super().save_model(request, obj, form, change)
            
            moleculesInfos = DataPackageManager.ExtractData(obj)
            
            try:
                for parentFile in moleculesInfos:
                    sortedMolecules = sorted(parentFile.get("molecules"), key=attrgetter("h")) # if you do not sort it, groupby is no working properly
                    reactionStepGroups = groupby(sortedMolecules, key=attrgetter("h"))
                    
                    for group in reactionStepGroups:
                        moleculesList = [] # This is for second method.                    
                        for molecules in group[1]:
                            moleculesList.append(molecules) 
    
                        # --- For determining parent molecules and eliminating unnecessary parent files  -----                    
                        parent = None
                        childMolecules = []
                        sortedMoleculesList = sorted(moleculesList, key=attrgetter("inchiKey")) # if you do not sort it, groupby is no working properly
                        inchiKeyGroups = groupby(sortedMoleculesList, key=attrgetter("inchiKey"))
                 
                        for mols in inchiKeyGroups:
                            ms = list(mols[1])
                            molsLen = len(ms)
                            if (molsLen != 1 and parent == None):
                                ms[0].isParent = True
                                ms[0].funcGroup = ""
                                if hasattr(ms[0], 'OptimizationResults'):
                                    if ms[0].OptimizationResults != None:
                                        for item in ms[0].OptimizationResults:
                                            item["isParent"] = True
                                            item["funcGroup"] = ""
                                if hasattr(ms[0], 'SingleResults'):
                                    if ms[0].SingleResults != None:
                                        for item in ms[0].SingleResults:
                                            item["isParent"] = True
                                            item["funcGroup"] = ""
                                parent = ms[0]
                            else:
                                ms[0].isParent = False
                                childMolecules.append(ms[0])
                            
                        try:
                            parentMolecule, created = Molecule.objects.get_or_create(inchiKey = parent.inchiKey, defaults={'smiles': parent.smiles})
                            parentMolecule, created = Molecule.objects.update_or_create(id = parentMolecule.id, defaults={'parentMolecule': parentMolecule})
                            DataPackageAdmin.save_calc(parent, parentMolecule, obj)
                        except Exception as e:
                            print("******************************************************************")
                            print("An error accured at database saving part. (parent)")
                            print(parent.filePath)
                            print(parent.actualName)
                            print(parent.__dict__)
                            if hasattr(e, 'message'):
                                print(e.message)
                            else:
                                print(e)
                            print("******************************************************************")
                            
                        for child in childMolecules: 
                            try:
                                childMolecule, createdChild = Molecule.objects.get_or_create(inchiKey = child.inchiKey, defaults={'smiles': child.smiles, 'parentMolecule': parentMolecule})
                                DataPackageAdmin.save_calc(child, childMolecule, obj)
                            except Exception as e:
                                print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
                                print("An error accured at database saving part. (parent)")
                                print(parent.filePath)
                                print(parent.actualName)
                                print(parent.__dict__)
                                print("An error accured at database saving part. (child)")
                                print(child.filePath)
                                print(child.actualName)
                                print(child.__dict__)
                                if hasattr(e, 'message'):
                                    print(e.message)
                                else:
                                    print(e) 
                                print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
            except:
                print("An errror accured after datapackage.")  
        else:
            print("DO SOMETHING")
     

class PairMatchAdmin(admin.ModelAdmin):
    def save_model(self, request, obj, form, change):
        super().save_model(request, obj, form, change)
        
        results = PairMatchManager.ExtractData(obj)    
        
        sortedPairList = sorted(results, key=attrgetter("inchiKey"))
        groupedPairList = groupby(sortedPairList, key=attrgetter("inchiKey"))
        parentPairInchiKeyList = []
        for p in groupedPairList:
            pList = list(p[1])
            if len(pList) != 1:
                parentPairInchiKeyList.append(p[0])
        
        for parentPair in parentPairInchiKeyList:
            for pair in results:
                if pair.inchiKey == parentPair:
                    pair.funcGroup =""
                          
        for pair in results:
            reactantMolecule = None
            productMolecule = None
            reactantSolutionEnergy = None
            productSolutionEnergy = None
            try:
                fucGroup = FunctionalGroup.objects.get(stoichiometry= pair.funcGroup)
            except FunctionalGroup.DoesNotExist:
                fucGroup = None
            
            try:
                productMolecule = Molecule.objects.get(inchiKey= pair.inchiKey)
            except Molecule.DoesNotExist:
                productMolecule = None
                productSolutionEnergy = None
                print("--- PRODUCT CAN NOT BE FOUND -------------------------------------------------------")
                print(pair.__dict__)
                print("------------------------------------------------------------------------------------")
            
            if productMolecule == None:
                splitedInchiKey = pair.inchiKey.split('-')[0]
                try:
                    productMolecule = Molecule.objects.get(inchiKey__startswith= splitedInchiKey)
                except Molecule.DoesNotExist:
                    productMolecule = None
                    productSolutionEnergy = None
                    
            try:
                productJob = Job.objects.filter(molecule= productMolecule, jobType='SinglePoint', dataPackage_id=obj.notes)
                if len(productJob) == 1:
                    productJob= productJob[0]
                elif len(productJob) == 0:
                    productJob=  None
                    productSolutionEnergy = None
                else:
                    lenOfProducts = len(productJob)
                    print("HATA")
                    print(productMolecule.__dict__)
                    print(productJob.__dict__)
                    print(lenOfProducts)
                    productJob= productJob[lenOfProducts-1]

                
                try:
                    if productJob is not None :
                        productSolutionEnergy = ScfCalc.objects.get(job= productJob.id)
                    else:
                        productSolutionEnergy = None
                except ScfCalc.DoesNotExist:
                    productSolutionEnergy = None
            except Job.DoesNotExist:
                productJob = None
                productSolutionEnergy = None
                
            # print("!!!!!!!!!!!!!!")
            # print(pair.actualName)
            # print(obj.notes)
            # print(fucGroup.id)
            # print(int(pair.h)-1)
            # print("!!!!!!!!!!!!!!")
            try:
                job = Job.objects.filter(name__startswith= pair.actualName, dataPackage_id=obj.notes, jobType = 'SinglePoint', functionalGroup_id = fucGroup.id, reactionStep = int(pair.h)-1)
                for item in job:
                    try:
                        reactantMolecule = Molecule.objects.get(id=item.molecule_id)
                        try:
                            reactantSolutionEnergy = ScfCalc.objects.get(job= item)
                        except ScfCalc.DoesNotExist:
                            reactantSolutionEnergy = None
                        break
                    except Molecule.DoesNotExist:
                        reactantMolecule = None
                        reactantSolutionEnergy = None
            except Job.DoesNotExist:
                job = None
                reactantMolecule = None
                reactantSolutionEnergy = None
            
            # print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            # print(reactantMolecule.__dict__)
            # print("__________________________________________________")
            # print(productMolecule.__dict__)
            # print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            
            # print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
            # print(reactantSolutionEnergy.__dict__)
            # print(productSolutionEnergy.__dict__)
            # print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
            
            
                
            
            
            if reactantMolecule != None and productMolecule != None:
                # print("++")
                reaction, reacCreated = Reaction.objects.get_or_create(
                    reactant = reactantMolecule,
                    product = productMolecule,
                    bondType = pair.reac,
                    pairPackage = obj)
                if reacCreated:
                    if reactantSolutionEnergy != None and productSolutionEnergy != None:
                        updatedReaction = Reaction.objects.filter(id = reaction.id).update(productEnergy=Decimal(productSolutionEnergy.solutionEnergy), reactantEnergy=Decimal(reactantSolutionEnergy.solutionEnergy), reactionEnergy=(Decimal(productSolutionEnergy.solutionEnergy)-Decimal(reactantSolutionEnergy.solutionEnergy)-Decimal(-1.16273669344)))
                    else:
                        if reactantSolutionEnergy != None:
                            updatedReaction = Reaction.objects.filter(id = reaction.id).update(reactantEnergy=reactantSolutionEnergy.solutionEnergy)
                        if productSolutionEnergy != None:
                            updatedReaction = Reaction.objects.filter(id = reaction.id).update(productEnergy=productSolutionEnergy.solutionEnergy)


                            
admin.site.register(DataPackage, DataPackageAdmin)
admin.site.register(PairMatchSmilesPackage, PairMatchAdmin)

