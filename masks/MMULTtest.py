import unittest
from MMULT import MMULT
class MMULTTEst(unittest.TestCase):
    def test1Load(self):
        mmult = MMULT("../../alignmentFiles/mmult.list-FINAL.aln")
        mmult.load()
        keys = mmult.alignments.keys()
        keys.sort()
        for key in keys:
            print("%s %s \n" % (key, mmult.alignments[key]))
            
        mmult.serialize("../../alignmentFiles/mmult.list-FINAL.alnF")
        
    def test2computePercentageOfAlignment(self):
        mmult = MMULT("../../alignmentFiles/mmult.list-FINAL.aln")
        mmult.load()
        mmult.computePercentageOfAlignment()
        for residue in mmult.residues.items():
            print ("[%d] = %f\n" % (residue[0], residue[1]))
        
    def test4findFirstAlignedResidue(self):
        mmult = MMULT("../../alignmentFiles/mmult.list-FINAL.aln")
        mmult.load()
        assert mmult.findFirstAlignedResidue("1a0hA02.pdb") == 0
        assert mmult.findFirstAlignedResidue("1avgL00.pdb") == 29
    
    def test5loadResidues(self):
        mmult = MMULT("../../alignmentFiles/mmult.list-FINAL.aln")
        mmult.load()
        mmult.loadResidues()
        residues = mmult.pdbResiduesOrder
        for (pdbFileName, number) in residues.items():
            print("%s - %d \n" %(pdbFileName, number))
    def test6getResidueNumber(self):
        mmult = MMULT("../../alignmentFiles/mmult.aux")
        mmult.load()
        pdbFileName  = "1a0hA02.pdb"
        pdbFileName2 = "1avgL00.pdb"
        pdbFileName3 = "1bbrJ00.pdb"
        assert  0 == mmult.getPositionInPDB(pdbFileName, 0)
        assert -1 == mmult.getPositionInPDB(pdbFileName, 1)
        assert  1 == mmult.getPositionInPDB(pdbFileName, 2)
        assert  2 == mmult.getPositionInPDB(pdbFileName, 3)
        assert  3 == mmult.getPositionInPDB(pdbFileName, 4)
        
        assert -1 == mmult.getPositionInPDB(pdbFileName2, 0)
        assert -1 == mmult.getPositionInPDB(pdbFileName2, 1)
        assert  0 == mmult.getPositionInPDB(pdbFileName2, 2)
        assert  1 == mmult.getPositionInPDB(pdbFileName2, 3)
        assert  2 == mmult.getPositionInPDB(pdbFileName2, 4)
        
        assert -1 == mmult.getPositionInPDB(pdbFileName3, 0)
        assert -1 == mmult.getPositionInPDB(pdbFileName3, 1)
        assert  0 == mmult.getPositionInPDB(pdbFileName2, 2)
        assert  1 == mmult.getPositionInPDB(pdbFileName2, 3)
        assert  2 == mmult.getPositionInPDB(pdbFileName2, 4)
        assert -1 == mmult.getPositionInPDB(pdbFileName3, 5)
        assert  3 == mmult.getPositionInPDB(pdbFileName3, 6)
        mmult = MMULT("../../alignmentFiles/mmult.list-FINAL.aln")
        mmult.load()
        assert 69  == mmult.getPositionInPDB(pdbFileName, 69)
        assert -1  == mmult.getPositionInPDB(pdbFileName2, 69)
        assert 34  == mmult.getPositionInPDB(pdbFileName3, 68)
        
        #print("%d"%(mmult.getResidueNumber(pdbFileName2, 2),))
        #assert  0 == mmult.getResidueNumber(pdbFileName2, 2)
        #assert  1 == mmult.getResidueNumber(pdbFileName2, 3)
        #assert -1 == mmult.getResidueNumber(pdbFileName3, 5)
    def test7percentageFileName(self):
        mmult = MMULT("../../alignmentFiles/mmult.list-FINAL.aln")
        percentageFileName = "../../alignmentFiles/percentage.aux"
        mmult.load()
        mmult.loadResidues()
        mmult.computePercentageOfAlignment()
        mmult.serializePercentageOfAlignment(percentageFileName)
