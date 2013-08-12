from entities.atom import Atom
from decimal import *
from numpy import matrix
import unittest

class AtomTestCase(unittest.TestCase):
    rawAtom  = "ATOM      1  N   SER A   4      12.140  45.657  29.217  1.00 40.81           N  "
    rawAtom2 = "ATOM   1162  CG2 VAL B 147X     10.793   0.886  28.197  1.00 11.51           C2+"
    rawAtom3 = "ATOM     48  N   PRO C   8      -2.513 -13.594  63.478  1.00  2.00           N"  
    rawAtom4 = "ATOM     49  CA  PRO D   8      -1.194 -14.125  63.727  1.00  2.00           C"  
    rawAtom5 = "ATOM     50  C   PRO E   8      -0.065 -13.490  62.933  1.00  6.90           C"  
    rawAtom6 = "ATOM     51  O   PRO F   8      -0.289 -12.746  61.987  1.00 16.83           O"  
    rawAtom7 = "ATOM     52  CB  PRO G   8      -1.402 -15.549  63.449  1.00  2.00           C"  
    def setUp(self):
        self.atom  = Atom(AtomTestCase.rawAtom)
        self.atom2 = Atom(AtomTestCase.rawAtom2) 
        self.atom3 = Atom(AtomTestCase.rawAtom3)
        self.atom4 = Atom(AtomTestCase.rawAtom4)
        self.atom5 = Atom(AtomTestCase.rawAtom5)
        self.atom6 = Atom(AtomTestCase.rawAtom6)
        self.atom7 = Atom(AtomTestCase.rawAtom7)
    def tearDown(self):
        self.atom  = None
        self.atom2 = None
    def test1ToString(self):
        assert AtomTestCase.rawAtom   == self.atom.toString()
        assert AtomTestCase.rawAtom2  == self.atom2.toString()
    def test2Shift(self):
        shift = matrix([0.01, 0.01, 0.01]).T
        self.atom.shift(shift)
        assert "ATOM      1  N   SER A   4      12.150  45.667  29.227  1.00 40.81           N  " == self.atom.toString()
    def test3Center(self):
        center = matrix([0.00, 0.00, 0.00]).T
        self.atom.center(center)
        assert AtomTestCase.rawAtom   == self.atom.toString()
    def test4inBackbone(self):
        assert self.atom.inBackbone()  == True
        assert self.atom3.inBackbone() == True
        assert self.atom4.inBackbone() == True
        assert self.atom6.inBackbone() == True
        assert self.atom7.inBackbone() == True
        
        assert self.atom2.inBackbone() == False
        assert self.atom5.inBackbone() == False
    
    def test5InChain(self):
        assert self.atom.inChain(" A ")
        assert self.atom2.inChain(" A ") == False
    
    def test6InSegment(self):
        assert self.atom.inSegment(4, 5)
        assert self.atom2.inSegment(148, 150) == False
        assert self.atom2.inSegment(147, 150) == True
    
    def test7atomWeight(self):
        assert self.atom.atomicWeight() == 14.0067
        
    def testinInsertionsSegment(self):
        rawStartAtom        = "ATOM      1  N   ALA L   1T     95.966  23.494 131.088  1.00 78.27           N  "
        rawInResidueAtom    = "ATOM     21  C   GLU L   1Q     90.085  26.425 128.760  1.00 63.99           C  "
        rawMiddleAtom       = "ATOM    158  N   CYS L   1      81.822  27.897 122.945  1.00 35.25           N  "
        rawEndAtom          = "ATOM    271  N   LYS L  14A     86.041  46.599 116.234  1.00 45.79           N  " 
        rawoutOfsegmentAtom = "ATOM    378  O   GLY L  14N     92.287  43.790 140.820  1.00101.24           O  " 
        rawoutOfsegmentAtom2= "ATOM    378  O   GLY L  15N     92.287  43.790 140.820  1.00101.24           O  "
        startAtom = Atom(rawStartAtom)
        inResidueAtom = Atom(rawInResidueAtom)
        middleAtom = Atom(rawMiddleAtom)
        endAtom = Atom(rawEndAtom)
        outOfsegmentAtom = Atom(rawoutOfsegmentAtom) 
        outOfsegmentAtom2 = Atom(rawoutOfsegmentAtom2)
        start = 1
        stop  = 14
        startInsertion = 'T'
        stopInsertion = 'M'
        startOrder = False
        stopOrder = True
        
        assert startAtom.inInsertionsSegment(start, stop, startInsertion, stopInsertion, startOrder, stopOrder) ==  True
        assert inResidueAtom.inInsertionsSegment(start, stop, startInsertion, stopInsertion, startOrder, stopOrder) ==  True
        assert middleAtom.inInsertionsSegment(start, stop, startInsertion, stopInsertion, startOrder, stopOrder) == True
        assert endAtom.inInsertionsSegment(start, stop, startInsertion, stopInsertion, startOrder, stopOrder) == True
        assert outOfsegmentAtom.inInsertionsSegment(start, stop, startInsertion, stopInsertion, startOrder, stopOrder) == False
        assert outOfsegmentAtom2.inInsertionsSegment(start, stop, startInsertion, stopInsertion, startOrder, stopOrder) == False
    def testgetResidue(self):
        assert self.atom.getResidue() == (4, ' ')
        assert self.atom2.getResidue() == (147,'X')