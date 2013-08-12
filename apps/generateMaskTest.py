import unittest
from apps.generateMask import *
class PDBTestCase(unittest.TestCase):
    def testOperate(self):
        workingDirectory  = "/home/jatienza/Desktop/cathAnalysis/eclipseProject/cathAnalysis/src/launchers/kk/"
        alnFileName       = workingDirectory+"1.10.10.180.mmult-FINAL.aln"
        minimumPercentage = 0.2
        verbose           = True
        pdbFileName       = "1aa7A02.pdb"
        resultPDBFileName = "1aa7A02_Results.pdb"
        operate(workingDirectory, alnFileName, minimumPercentage, pdbFileName, resultPDBFileName,verbose)