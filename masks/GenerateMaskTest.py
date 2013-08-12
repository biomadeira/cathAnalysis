import unittest
import GenerateMask
import os
class GenerateMaskTestCase(unittest.TestCase):
    def atest1(self):
        workingDirectory  = "../../alignmentFiles/"
        #alnFileName       = workingDirectory+"mmult3.aux"
        alnFileName       = workingDirectory+"mmult.list-FINAL.aln"
        #maskedPDBfileName = "1a0hA02.pdb"
        maskedPDBfileName = "1avgL00.pdb"
        loader         = GenerateMask.MaskDataLoader(workingDirectory,alnFileName) 
        loader.load(0.3)
        
        #residues = loader.getPDBresidues()
        #assert len(residues) > 0
        
        #print("len(residues) = %d\n" %(len(residues),))
        
        #for key, value in loader.getPDBresidues().items():            
        #    print("%d %d \n" %(key, value))
        
        loader.compute_masked_residues(maskedPDBfileName,True)
        resultFileName = workingDirectory+"mask.pdb"
        loader.serialize_masked_atoms(resultFileName)
        bildFileName = workingDirectory+"mask.bild"
        loader.serialize_masked_atoms(bildFileName, True)
        histogramFileName = workingDirectory + "histogram"
        loader.generateHistogramOfPercentageOfAlignment(histogramFileName,show=False,serializePercentageOfAlignment=True)
        loader.serialize_variance_results()
    def test2MaskGenerator(self):
        workingDirectory  = "../../alignmentFiles"
        alnFileName       = "mmult.list-FINAL.aln"
        varianceFileName  = "variance.results"
        maskFileName      = "varianceMask.results"
        pdbFileName       = "1avgL00.pdb"
        resultPDBFileName = "1avgL00_mask.pdb"
    
        maskGenerator = GenerateMask.MaskGenerator(workingDirectory, alnFileName, 0.3, varianceFileName, maskFileName, pdbFileName, resultPDBFileName)
        maskGenerator.load()
        maskGenerator.compute_masked_residues(False)
        maskGenerator.serialize_masked_atoms(True)