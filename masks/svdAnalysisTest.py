import unittest
import numpy
from  masks import svdAnalysis 
from masks import generateMultipleAlignment

class SVDAnalysisTestCase(unittest.TestCase):
    def test1loadCAmatrix(self):
        """To generate filteredResults2.Atomic:
        1.-generateMultipleAlignment.
            alnFileName = "../../alignmentFiles/mmult3.aux"
            resultsFileName2 = "../../alignmentFiles/results2.results"
            generateMultipleAlignment.getResidues(0.8, resultsFileName2, alnFileName, workingDirectory, False,True, False)
            AtomicresultsFileName = "../../alignmentFiles/results2.Atomic"
            generateMultipleAlignment.composeAtomicMatrix(resultsFileName2, AtomicresultsFileName)
            filteredAtomicMatrixFileName = "../../alignmentFiles/filteredResults2.Atomic"
            generateMultipleAlignment.filterAtomicMatrix(AtomicresultsFileName, filteredAtomicMatrixFileName)
         """
        caMatrixFileName = "../../alignmentFiles/filteredResults2.Atomic" 
        dimensions = generateMultipleAlignment.getDimensions(caMatrixFileName)
        (mask, cAmatrix) = svdAnalysis.loadCAmatrix(caMatrixFileName, dimensions)
        
        #assert (70,9) == mask.shape
        #assert (70,9) == cAmatrix.shape
        print ("\n")
        print cAmatrix
        print ("\n")
        print mask
        print ("\n")
        

    def test4matrixMean(self):
        mask = numpy.array([[True, True, True, True, True, True],[True, True, True, True, True, True],[True, True, True, True,True, True]])
        #mask = numpy.array([[False,False, False, True, True, True],[True, True, True, True, True, True],[True, True, True, True, True, True]])
        matrix = numpy.array([[0,1,2,3,4,5],[10,11,12,13,14,15],[20,21,22,23,24,25]], dtype=float)
        print matrix
        svdAnalysis.matrixMean(mask, matrix)
        print matrix

    def test4rowMeanAtom(self):
        mask   = numpy.array([True, True, True, True, True, True])
        matrix = numpy.array([0,1,2,3,4,5])

        mean = svdAnalysis.rowMeanAtom(mask, matrix)
        
        print mean
        
    def test5getCovarianceMatrix(self):
        mask = numpy.array([[False,False, False, True, True, True],[True, True, True, True, True, True]])
        #mask = numpy.array([[True, True, True, True, True, True],[True, True, True, True, True, True]])
        matrix = numpy.array([[0,1,2,3,4,5],[10,11,12,13,14,15]], dtype=float)
        (caMatrix,caMask) = svdAnalysis.getCovarianceMatrix(mask, matrix, True)
        print("\n caMatrix\n")
        print caMatrix
        print("\ncaMask\n")
        print caMask
        (ca_atoms, ca_mask) = svdAnalysis.getAtomsVariance(caMatrix, caMask)
        
    def testGlobal(self):
        caMatrixFileName = "../../alignmentFiles/filteredResults2.Atomic"
        dimensions = generateMultipleAlignment.getDimensions(caMatrixFileName)
        print ("%d %d\n"%(dimensions[0],dimensions[1]))
        (mask, matrix)   = svdAnalysis.loadCAmatrix(caMatrixFileName, dimensions)  
        (caMatrix,caMask) = svdAnalysis.getCovarianceMatrix(mask, matrix)
        print ("%d %d\n"%(caMatrix.shape[0],caMatrix.shape[1]))
        matrixFileName = "../../alignmentFiles/matrix.float" 
        maskFileName  = "../../alignmentFiles/mask.float"
        #print("\n caMatrix\n")
        #print caMatrix
        #print("\ncaMask\n")
        #print caMask
        svdAnalysis.serialize_matrix(caMatrix, matrixFileName)
        diagonalFileName = "../../alignmentFiles/diagonal.float"
        (ca_atoms, ca_mask) = svdAnalysis.getAtomsVariance(caMatrix, caMask)
        svdAnalysis.serialize_matrix(ca_atoms, diagonalFileName)
        svdAnalysis.serialize_matrix(ca_mask, maskFileName)
        
    
    def testGlobalReal(self):
        alnFileName = "../../alignmentFiles/mmult3.aux"
        workingDirectory = "../../alignmentFiles/"
        resultsFileName2 = "../../alignmentFiles/results2.results"
        generateMultipleAlignment.getResidues(0.4, resultsFileName2, alnFileName, workingDirectory, False,True, False)
        AtomicresultsFileName = "../../alignmentFiles/results2.Atomic"
        generateMultipleAlignment.composeAtomicMatrix(resultsFileName2, AtomicresultsFileName)
        filteredAtomicMatrixFileName = "../../alignmentFiles/filteredResults2.Atomic"
        generateMultipleAlignment.filterAtomicMatrix(AtomicresultsFileName, filteredAtomicMatrixFileName)
        caMatrixFileName = "../../alignmentFiles/filteredResults2.Atomic"
        dimensions = generateMultipleAlignment.getDimensions(caMatrixFileName)
        print "dimensions\n"
        print ("%d %d \n"%(dimensions[0],dimensions[1]))
        (mask, matrix)   = svdAnalysis.loadCAmatrix(caMatrixFileName, dimensions)  
        (caMatrix,caMask) = svdAnalysis.getCovarianceMatrix(mask, matrix)
        print ("matrix.shape %d %d" % (caMatrix.shape[0],caMatrix.shape[1]))
        maskFileName   = workingDirectory +  "mask.atomic"
        matrixFileName = workingDirectory + "matrix.atomic"
        svdAnalysis.serializeCAMaskAndCAMatrix(mask, matrix, maskFileName, matrixFileName)
        realCAMatrixFileName = workingDirectory + "caMatrix.atomic"
        realCAMaskFileName   = workingDirectory + "caMask.atomic"
        svdAnalysis.serializeCAMaskAndCAMatrix(caMask, caMatrix, realCAMaskFileName, realCAMatrixFileName)
        
        
        
        
    def testStudyVariance(self):
        alnFileName = "../../alignmentFiles/mmult.list-FINAL.aln"
        workingDirectory = "../../alignmentFiles/"
        resultsFileName2 = "../../alignmentFiles/results2.results"
        generateMultipleAlignment.getResidues(0.8, resultsFileName2, alnFileName, workingDirectory, False,True, False)
        AtomicresultsFileName = "../../alignmentFiles/results2.Atomic"
        generateMultipleAlignment.composeAtomicMatrix(resultsFileName2, AtomicresultsFileName)
        filteredAtomicMatrixFileName = "../../alignmentFiles/filteredResults2.Atomic"
        generateMultipleAlignment.filterAtomicMatrix(AtomicresultsFileName, filteredAtomicMatrixFileName)
        caMatrixFileName = "../../alignmentFiles/filteredResults2.Atomic"
        dimensions = generateMultipleAlignment.getDimensions(caMatrixFileName)
        (mask, matrix)   = svdAnalysis.loadCAmatrix(caMatrixFileName, dimensions)  
        (caMatrix,caMask) = svdAnalysis.getCovarianceMatrix(mask, matrix)
        (atoms_variance,atoms_variance_mask) = svdAnalysis.getAtomsVariance(caMatrix,caMask)
        (max, min, mean, normalizedAtoms) = svdAnalysis.atoms_variance_statistics(atoms_variance,atoms_variance_mask)
        print ("max %f \n" %(max,))
        print ("min %f \n" %(min,))
        print ("mean %f \n" %(mean,))
        print normalizedAtoms
        
        print "\n\n"
        workingDirectory = "../../alignmentFiles/"
        varianceFileName = workingDirectory+"variance.atomic"
        maskFileName     = workingDirectory+"variance_mask.atomic"
        svdAnalysis.serializeAtomsVariance(atoms_variance, atoms_variance_mask, varianceFileName, maskFileName)
        (atoms_variance,atoms_variance_mask) = svdAnalysis.encarnateAtomsVariance(varianceFileName, maskFileName)
        assert 30 == len(atoms_variance[0])
        for variance in atoms_variance[0]:
            print ("%3.f \n"%(variance,))
        
        
        
        
                
        
        
        