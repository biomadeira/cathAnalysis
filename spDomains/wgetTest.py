import unittest
import os
import spDomains.wget
class PDBTestCase(unittest.TestCase):
    def atest1compress(self):
        workingDirectory = "/home/jatienza/.eclipse/workspace/cath/pdbs"
        fileNames =[]
        fileNames.append(os.path.join(workingDirectory,"1a0hA02.pdb"))
        fileNames.append(os.path.join(workingDirectory,"1a0hD02.pdb"))
        compressedFileName = os.path.join(workingDirectory, "tt")
        compressFile(fileNames, compressedFileName)
        
    def atest2INB(self):
        pdbCode = "1a0h"
        directory="/home/inab/databases/flat/pdb"
        fileName = "/home/inab/databases/flat/pdb/pdb1a0h.ent.Z"
        
        assert fileName == spDomains.wget.wgetINBPDBFileName(pdbCode,directory)
        
    def testUnzipINBPDB(self):
        inputFileName  = "/home/jatienza/.eclipse/workspace/cath/pdbs/pdb1tbr.ent.Z"
        outputFileName = "/home/jatienza/.eclipse/workspace/cath/pdbs/1tbr_inb.pdb"
        spDomains.wget.uncompressZCATfiles(inputFileName, outputFileName)
        
        