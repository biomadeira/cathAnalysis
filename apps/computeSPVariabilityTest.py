import unittest
from apps.computeSPVariability2 import *
class PDBTestCase(unittest.TestCase):
    def testOperate(self):
        workingDirectory  = "/home/jatienza/Desktop/cathAnalysis/eclipseProject/cathAnalysis/src/launchers/kk/"
        alnFileName       = workingDirectory+"1.10.10.180.mmult-FINAL.aln"
        minimumPercentage = 0.2
        verbose           = True
        operate(workingDirectory, alnFileName, minimumPercentage, verbose)