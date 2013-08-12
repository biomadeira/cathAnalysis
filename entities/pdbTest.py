from entities.pdb import PDB
from entities.atom import Atom
from numpy import matrix
from libraries import rotationMatrices
from math import radians
from string import rjust, join

import unittest
class PDBTestCase(unittest.TestCase):
    pdbFileName = "../../pdbs/1CBO.pdb"
    
    def test01RASH(self):
        queryFileName = "../../pdbs/1cf3A03.pdb"
        pdb = PDB(queryFileName)
        pdb.load_atoms_hetams()
        center = matrix([42.639, 21.304, 46.057]).T
        shift  = matrix([15.155, -0.310, -36.823]).T
        rotationMatrix = matrix ([[-0.030, 0.982, 0.186], [-0.999, -0.036, 0.029], [0.035, -0.185, 0.982]])
        pdb.ash_transform(center, rotationMatrix, shift)
        
        pdb.serialize("../../pdbs/1cf3A03_RT.pdb")
        
        templateFileName = "../../pdbs/1gpeA03.pdb"
        pdb = PDB(templateFileName)
        pdb.load_atoms_hetams()
        center = matrix([42.639, 21.304, 46.057]).T
        shift  = matrix([15.155, -0.310, -36.823]).T
        rotationMatrix = matrix ([[-0.030, 0.982, 0.186], [-0.999, -0.036, 0.029], [0.035, -0.185, 0.982]])
        pdb.ash_transform(center, rotationMatrix, shift)
        
        pdb.serialize("../../pdbs/1gpeA03_RT.pdb")
    
    def test02BasicOperations(self):
        pdbFileName = "../../pdbs/1.pdb"
        pdb = PDB(pdbFileName)
        pdb.load_atoms_hetams()
        x_shift  = matrix([10, 0, 0]).T
        y_shift  = matrix([0, 10, 0]).T
        z_shift  = matrix([0, 0, 10]).T
        
        pdb.shift(x_shift);
        pdb.serialize("../../pdbs/x_shift.pdb")
        
        pdb = PDB(pdbFileName)
        pdb.load_atoms_hetams()
        pdb.shift(y_shift)
        pdb.serialize("../../pdbs/y_shift.pdb")
        
        pdb = PDB(pdbFileName)
        pdb.load_atoms_hetams()
        pdb.shift(z_shift)
        pdb.serialize("../../pdbs/z_shift.pdb")
        
        
        alpha =  radians(45)
        rotationMatrix = rotationMatrices.XYZ(alpha, alpha, alpha)
        center = matrix([0, 0, 0]).T
        pdb = PDB(pdbFileName)
        pdb.load_atoms_hetams()
        pdb.rotate(rotationMatrix, center)
        pdb.serialize("../../pdbs/R_xyz45_45_45.pdb")
        
        alpha =  radians(45)
        rotationMatrix = rotationMatrices.ZXY(alpha, alpha, alpha)
        center = matrix([0, 0, 0]).T
        pdb = PDB(pdbFileName)
        pdb.load_atoms_hetams()
        pdb.rotate(rotationMatrix, center)
        pdb.serialize("../../pdbs/R_zxy45_45_45.pdb")
        
        alpha =  radians(45)
        rotationMatrix = rotationMatrices.YXZ(alpha, alpha, alpha)
        center = matrix([0, 0, 0]).T
        pdb = PDB(pdbFileName)
        pdb.load_atoms_hetams()
        pdb.rotate(rotationMatrix, center)
        pdb.serialize("../../pdbs/R_yxz45_45_45.pdb")
        
        
    def test03Transform(self):
        self.pdb = PDB(PDBTestCase.pdbFileName)
        self.pdb.load_atoms_hetams()
        center = matrix([0, 0, 0]).T
        shift  = matrix([1, 1, 1]).T
        rotationMatrix = rotationMatrices.XYZ(radians(90), radians(90), radians(90))
        self.pdb.chimera_transform(center, rotationMatrix, shift)
        self.pdb.serialize("../../pdbs/transformed.pdb")
        self.pdb = None   
    def test04Itk_ops(self):
        self.pdb = PDB(PDBTestCase.pdbFileName)
        self.pdb.load_atoms_hetams()
        rotationMatrix = matrix ([[0.742755, -0.0816224, 0.664569], [0.286523, 0.93582, -0.205295], [-0.605161, 0.342898, 0.718472]])
        shift = matrix([-6.01415, 4.92925, 16.1776]).T
        """shift = matrix([2,5,4]).T"""
        target_first_voxel = matrix([58, -14, 28]).T    
        model_first_voxel = matrix([60, -14, 28]).T
        self.pdb.itk_transform(rotationMatrix, shift)
        print "completed"
        self.pdb.serialize("../../pdbs/1jlvC1.pdb")
        self.pdb = None
    
    def test05getAtomsFromSegment(self):
        pdbFileName = "../../pdbs/1A00.pdb"
        chain = " A "
        start = 1
        stop  = 141
        startInsertion = ''
        stopInsertion = ''
        startOrder = True
        stopOrder = True
        pdb = PDB(pdbFileName)
        pdb.load_atoms_hetams()
        
        segments = pdb.getAtomsFromSegment(chain, start, stop, startInsertion, stopInsertion, startOrder, stopOrder)
        assert segments is not None
        
        for (lineNumber, atom) in segments.items():
            print "[%d] = %s" %(lineNumber, atom.toString())
    
    def test09getAtomsFromSegment2(self):
        pdbFileName = "../../pdbs/2J9N.pdb"
        segmentFileName = "../../pdbs/2J9N_A_101_102.pdb"
        chain = " A "
        start = 101
        stop  = 102
        startInsertion = ''
        stopInsertion = ''
        startOrder = True
        stopOrder = True
        pdb = PDB(pdbFileName)
        pdb.load_atoms_hetams()
        segment =[]
        #chain, start, stop, startInsertion, stopInsertion, startOder, stopOrder, model=0, onlyBackbone=False):
        segment = pdb.getAtomsFromSegment(chain, start, stop, startInsertion, stopInsertion, startOrder, stopOrder)
        assert segment is not None
        
        file = open(segmentFileName,"w")
        orderedKeys = segment.keys()
        orderedKeys.sort()
        
        for lineNumber in orderedKeys:
                file.write("%s\n" %(segment[lineNumber].toString(),))
        file.write("%s\n" % ("END"+77*" ",))
        file.close()
        
            
    def test06getAtomsFromDomain(self):
        pdbFileName = "../../pdbs/1CBO.pdb"
        domainFileName = "../../pdbs/1cboA01.pdb"
        file = open(domainFileName, "w")
        pdb = PDB(pdbFileName)
        pdb.load_atoms_hetams()
        segments =[]
        segments.append(('A', 9, 11,'',''))
        #segments.append(('A',385,401))
        #segments.append(('A',443,505))
        
        domains = pdb.getAtomsFromDomain(segments)
        
        for segment in domains:
            orderedKeys = segment.keys()
            orderedKeys.sort()
            for lineNumber in orderedKeys:
                file.write("%s\n" %(segment[lineNumber].toString(),))
        file.close()
    def test13getAtomsFromDomain(self):
        pdbFileName = "../../pdbs/1JMO.pdb"
        domainFileName = "../../pdbs/1jmoL00.pdb"
        file = open(domainFileName, "w")
        pdb = PDB(pdbFileName)
        pdb.load_atoms_hetams()
        segments =[]
        segments.append(('L', 1, 14,'B','M'))
        #segments.append(('A',385,401))
        #segments.append(('A',443,505))
        
        domains = pdb.getAtomsFromDomain(segments)
        
        for segment in domains:
            orderedKeys = segment.keys()
            orderedKeys.sort()
            for lineNumber in orderedKeys:
                file.write("%s\n" %(segment[lineNumber].toString(),))
        file.close()
    def test07CenterOfMass(self):
        pdbFileName = "../../pdbs/1cboA01.pdb"
        pdb = PDB(pdbFileName)
        pdb.load_atoms_hetams()
        print "center of mass: %f %f %f\n" % (pdb.centerOfMass()[0], pdb.centerOfMass()[1], pdb.centerOfMass()[2])
        
    def test08highestOccupancyLevel(self):
        rawAtom = "ATOM    661  CA AMET A 101      16.070  45.718  27.607  0.20  7.83           C "
        anAtom = Atom(rawAtom)
        pdbFileName = "../../pdbs/multipleLocations.pdb"
        pdb = PDB(pdbFileName)
        pdb.load_atoms_hetams()
        assert 'A' == anAtom.alternateLocation()
        print "%c" % (pdb.getHighestOccupancyLocation(anAtom),)
        assert 'B' == pdb.getHighestOccupancyLocation(anAtom)
        
    def test10orderOfInsertions(self):
        pdbFileName = "../../pdbs/1JMO.pdb"
        pdb = PDB(pdbFileName)
        pdb.load_atoms_hetams()
        assert pdb.orderOfInsertions(1) == False
        assert pdb.orderOfInsertions(14) == True
    
    def test14Fasta(self):
        pdbFileName = "../../pdbs/1JMO.pdb"
        pdb = PDB(pdbFileName)
        pdb.load_atoms_hetams()
        fastaFileName = "../../pdbs/1JMO.fasta"
        pdb.fasta(fastaFileName)
    
    def atest15Fasta2(self):

        pdb1avgL00 = "../../pdbs/1avgL00.pdb"
        pdb1bbrJ00 = "../../pdbs/1bbrJ00.pdb"
        fasta1avgL00 = "../../pdbs/1avgL00.fasta"
        fasta1bbrJ00 = "../../pdbs/1bbrJ00.fasta"
        
        Pdb1avgL00 = PDB(pdb1avgL00)
        Pdb1bbrJ00 = PDB(pdb1bbrJ00)
        
        Pdb1avgL00.load_atoms_hetams()
        Pdb1bbrJ00.load_atoms_hetams()
        
        Pdb1avgL00.fasta(fasta1avgL00)
        Pdb1bbrJ00.fasta(fasta1bbrJ00)
        
    def test16loadResidues(self):
        pdbFileName = "../../pdbs/1tbqL00.pdb"
        pdb = PDB(pdbFileName)
        pdb.load_atoms_hetams()
        pdb.loadResidues()
        assert len(pdb.residues) > 0
                              

        
    def test17getFirstResidue(self):
        
        pdblistFileName = "../../pdbs/pdbFileList.lst"
        workingDirectory = "../../pdbs/"
        pdListFile = open(pdblistFileName,"r")
        for pdbFileName in pdListFile:
            pdb = PDB(workingDirectory+pdbFileName.split()[0])
            pdb.load_atoms_hetams()
            firstResidue = pdb.getFirstResidue()
            print ("%s = (%d % c) \n" % (pdbFileName,firstResidue[0],firstResidue[1]))
        
    def test18getIntSeqRes(self):
        pdbFileName = "../../pdbs/1tbqL00.pdb"
        pdbFileName2 = "../../pdbs/1a0hA02.pdb"
        pdb = PDB(pdbFileName)
        pdb2 = PDB(pdbFileName2)
        pdb.load_atoms_hetams()
        pdb2.load_atoms_hetams()
        pdb.loadResidues()
        pdb2.loadResidues()
        pdb.loadIntSeqRes()
        pdb2.loadIntSeqRes()
        assert pdb.getIntSeqRes(0) == (1,'U')
        assert pdb.getIntSeqRes(1) == (1,'T')
        
        assert pdb2.getIntSeqRes(0) == (253,' ')
                              
    def test19getAtomsFromResidue(self):
        pdbFileName = "../../pdbs/1tbqL00.pdb"
        pdb  = PDB(pdbFileName)
        pdb.load_atoms_hetams()
        pdb.loadResidues()
        pdb.loadIntSeqRes()
        atoms = pdb.getAtomsFromResidue(0)
        assert len(atoms) > 0
        assert len(atoms) == 7
        for atom in atoms:
            print ("%s\n" % (atom.toString(),))
        
        #    print("%.3f %.3f %.3f" % (atom[0,0],atom[1,0], atom[2,0]))
        
    def test20sphereRadius(self):
        pdbFileName = "../../pdbs/1tbqL00.pdb"
        pdb  = PDB(pdbFileName)
        pdb.load_atoms_hetams()
        (cm, radius) = pdb.sphereRadius()
        print ("radius %f \n"%(radius,))
        
        sb = []
        sb.extend(rjust("%.3f" % cm[0],8))
        sb.extend(rjust("%.3f" % cm[1],8))
        sb.extend(rjust("%.3f" % cm[2],8))
        print ("%s\n", join(sb,''))
        
    def test21getCAlphaFromResidue(self):
        pdbFileName = "../../pdbs/1a0hA02.pdb"
        pdb  = PDB(pdbFileName)
        pdb.load_atoms_hetams()
        pdb.loadResidues()
        pdb.loadIntSeqRes()
        atom = pdb.getCAlphaFromResidue(0)
        assert atom.toString() == "ATOM      2  CA  ASP A 253      78.545  46.395   3.454  1.00 61.74           C  "
        atom = pdb.getCAlphaFromResidue(12)
        assert atom.toString() == "ATOM     89  CA  ASP A 264A     63.823  29.490 -10.700  1.00 86.06           C  "
