import unittest
from launcher import *
class LauncherTestCase(unittest.TestCase):
    def testOperate(self):
        launcher = Launcher("/home/jatienza/Desktop/cathAnalysis/eclipseProject/cathAnalysis/src/launchers/kk/sflexfit.lch")
        launcher.operate()
        
        