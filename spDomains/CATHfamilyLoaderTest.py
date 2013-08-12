import unittest
from spDomains.CATHfamilyLoader import CATHfamilyLoader
class CATHfamilyLoaderTest(unittest.TestCase):
    #cathcode = "2.60.40.670"
    #cathcode = "4.10.140.10"
    #cathcode = "3.30.560.10"
    #cathcode = "3.90.640.10"
    #cathcode = "4.10.140.10"
    #cathcode = "1.10.10.180"
    #host   = "okazaki"
    #user   = "sflexfit"
    #passwd = "M8z2LeGR"
    #db     = "sflexfit"
    cathcode = "1.10.10.410"
    host   = "localhost"
    user   = "sflexfit"
    passwd = "sflexfit"
    db     = "sflexfit"
    workingDirectory = "/home/jatienza/Desktop/cathExample/"
    CATHfamilyLoader = CATHfamilyLoader(cathcode,host,user, passwd,db, workingDirectory)
    
    def test1GetFamilyPDBS(self):
        self.CATHfamilyLoader.getFamilyPDBS()
        assert self.CATHfamilyLoader.pdbIDs is not None
    
    def test2downloadPDBfiles(self):
        self.CATHfamilyLoader.downloadPDBfiles()
    
    def test3UncompressFiles(self):
        self.CATHfamilyLoader.uncompressPDBFiles()
    
    def test4ExtractDomains(self):
        self.CATHfamilyLoader.extractDomains()
    
    def test5CompressResults(self):
        compressedFileName= "/home/jatienza/Desktop/cathExample"
        logFileName = self.CATHfamilyLoader.compressResults(compressedFileName)        
    def atest5Remove(self):
        self.CATHfamilyLoader.deleteTempFiles(verbose=True)
        self.CATHfamilyLoader.deleteAll(verbose=True)
    
    
