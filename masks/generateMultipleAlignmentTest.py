import unittest
import generateMultipleAlignment 
class GenerateMultipleAlignmentTestCase(unittest.TestCase):
    def test1getResidues(self):
        #alnFileName = "../../alignmentFiles/mmult.list-FINAL.aln"
        alnFileName = "../../alignmentFiles/mmult3.aux"
        resultsFileName  = "../../alignmentFiles/results.results"
        resultsFileName2 = "../../alignmentFiles/results2.results"
        workingDirectory = "../../alignmentFiles/"
        #generateMultipleAlignment.getResidues(0.8, resultsFileName, alnFileName, workingDirectory)
        generateMultipleAlignment.getResidues(0.4, resultsFileName2, alnFileName, workingDirectory, False,True, False)
        AtomicresultsFileName = "../../alignmentFiles/results2.Atomic"
        generateMultipleAlignment.composeAtomicMatrix(resultsFileName2, AtomicresultsFileName)
        filteredAtomicMatrixFileName = "../../alignmentFiles/filteredResults2.Atomic"
        generateMultipleAlignment.filterAtomicMatrix(AtomicresultsFileName, filteredAtomicMatrixFileName)
        
        finalAtomicMatrixFileName = "../../alignmentFiles/finalAtomic.Atomic"
        generateMultipleAlignment.getAtomicCoordinatesMatrix(filteredAtomicMatrixFileName, finalAtomicMatrixFileName)
        
        (rows, colums) = generateMultipleAlignment.getDimensions(filteredAtomicMatrixFileName)
        print ("rows: %d columns:%d\n" %(rows, colums))
        