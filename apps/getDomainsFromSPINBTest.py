import unittest
from apps.getDomainsFromSPINB import operate
class PDBTestCase(unittest.TestCase):
    def testOperate(self):
        inputDirectory = "/home/jatienza/tempWorkingDirectory/1.10.10.180"
        workingDirectory = "/home/jatienza/tempWorkingDirectory/1.10.10.180/output/"
        
        host   = "okazaki.cnb.csic.es"
        user   = "poli"
        passwd = "R7xvgAKK"
        db     = "sflexfit"
        cathcode = "1.10.10.180"
        operate(cathcode, inputDirectory, workingDirectory, host, user, passwd, db, verbose=True)