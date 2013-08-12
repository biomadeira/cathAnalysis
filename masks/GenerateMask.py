import os
import numpy
import math
import pylab
import MMULT

from masks import  generateMultipleAlignment

from masks import svdAnalysis

from entities import pdb

import libraries.OScommands

class MaskGenerator:
    colors = [".color blue", ".color dark green", ".color olive drab", ".color brown", ".color red", ".color orange red"]
    def __init__(self, workingDirectory, alnFileName, minimumPercentage, varianceFileName, maskFileName, pdbFileName, resultPDBFileName):

        self.__workingDirectory      = libraries.OScommands.gotoWorkingDirectory(True, workingDirectory)
        self.__alnFileName           = alnFileName
        self.__varianceFileName      = varianceFileName
        self.__maskFileName          = maskFileName
        self.__maskedPDBFileName     = pdbFileName
        self.__resultPDBFileName     = resultPDBFileName
        
        self.__minimumPercentage     = minimumPercentage
    def operate(self,verbose=False):
        """run the program with the default parameters """
        if verbose:
            print "loading data\n"
        self.load()
        if verbose:
            print "computing masked residues\n"
        self.compute_masked_residues(True)
        if verbose:
            print "serializing results\n"
        self.serialize_masked_atoms(False)
        
    def load(self):
        self.__load_atoms_variance()
        self.__load_mmult()
        
    def __load_mmult(self):
        self.__mmult = MMULT.MMULT(self.__alnFileName)
        self.__mmult.load()
        self.__mmult.computePercentageOfAlignment()
    
    def __load_atoms_variance(self):
        (self.atoms_variance, self.atoms_variance_mask) = svdAnalysis.encarnateAtomsVariance(self.__varianceFileName, self.__maskFileName)
 
    def compute_masked_residues(self, applyPercentage=False):
        """
        for each atom, give me the variance value of its residue in the multiple alignment
        stores in self.__masked_atoms a list with tuples (atom, variance,percentageOfAlignment) in order
        """
        self.__masked_atoms = []
        self.__maskedPDB = pdb.PDB(self.__maskedPDBFileName)
        self.__maskedPDB.load_atoms_hetams()
        self.__maskedPDB.loadResidues()
        self.__maskedPDB.loadIntSeqRes()
        
        self.__maskedPDBFileName = os.path.basename(self.__maskedPDBFileName)
        
        atom_order = 0
        for position_in_alignment in range(len(self.__mmult.residues)):
            if self.__mmult.residues[position_in_alignment] > self.__minimumPercentage:
                positionInPDB = self.__mmult.getPositionInPDB(self.__maskedPDBFileName, position_in_alignment)
                if positionInPDB != -1:
                    atom = self.__maskedPDB.getCAlphaFromResidue(positionInPDB)
                    if atom is not None:
                        variance = self.atoms_variance[0][atom_order]
                        if applyPercentage:
                           self.__masked_atoms.append((atom, self.__mmult.residues[position_in_alignment]*variance))
                        else:
                           self.__masked_atoms.append((atom, variance))
                atom_order+=1

    def serialize_masked_atoms(self, toBild=False):
        result = open(self.__resultPDBFileName, "w")
        if toBild:
            #(max, min, mean, normalizedAtoms) = svdAnalysis.atoms_variance_statistics(self.atoms_variance,self.atoms_variance_mask)
            self.computeColorsMap()
            for atom, variance in self.__masked_atoms:
                result.write(self.get_color(variance))
                result.write("\n")
                result.write(atom.toBildString())
                result.write("\n")
        else:
            for atom, variance in self.__masked_atoms:
                result.write(atom.maskedAtomToString(variance))
                result.write("\n")
                
        result.close()
        
    def computeColorsMap(self):
        differentValues = set()
        for atom, variance in self.__masked_atoms:
            differentValues.add(variance)
           
        groups = math.ceil(len(differentValues)/len(self.colors))
        l = []
        l.extend(differentValues)
        l.sort()
        i      = 0
        color  = 0
        self.__map = {}
        self.kk = self.__map
        for value in l:
            if i>groups:
               i=0
               try:
                   self.__map[value]=self.colors[color]
               except IndexError:
                    self.__map[value]=self.colors[len(self.colors)-1]
               color+=1
            else:
              i+=1
        
            
        self.__map[value]=self.colors[len(self.colors)-1]
        
        self.__maxKeys = self.__map.keys()
        self.__maxKeys.sort()

    def get_color(self, value):
        self.kkk = self.__maxKeys
        for i in range(len(self.__maxKeys)):
            if value < self.__maxKeys[i]:
                return self.__map[self.__maxKeys[i]] 
        
        return self.colors[len(self.colors)-1]

    


class MaskDataLoader:
    """1.-execute load with the desired percentage of variance
       2.-compute_masked_residues returns a list filled with sorted tuples (atom, variance) 
    """
    colors = [".color blue", ".color dark green", ".color olive drab", ".color brown", ".color red", ".color orange red"]
    def __init__(self, workingDirectory, alnFileName, patternPDBFileName='',bildFileName = '', histogramFileName=''): 
        self.__workingDirectory  = libraries.OScommands.gotoWorkingDirectory(True, workingDirectory)
        print ("wD:%s\n"%(self.__workingDirectory))
        self.__alnFileName       = os.path.join(self.__workingDirectory,alnFileName)
        self.__patternPDBFileName = patternPDBFileName
        self.__bildFileName = bildFileName
        self.__histogramFileName = histogramFileName
        
    def operate(self, minimumPercentage,stats=False, verbose=False):
        if verbose:
            print "loading data"
        self.load(minimumPercentage, verbose)
        if verbose:
            print "serializing variance results  \n"
        self.serialize_variance_results()
        
        if stats:
            if verbose:
                print "generating stats \n"
            self.compute_masked_residues(self.__patternPDBFileName)
            self.generate_stats(self.__bildFileName, self.__histogramFileName)
                
    def run_mmult_scripts(self):
        """@todo: must run the scripts to generate the alignment File (as returned from mmult)"""
        self.__alnFileName  = ""
        
    def __generate_file_names(self):
        """Generates the intermediate result filenames """
        self.__alphaCarbonsFileName          = "alphaCarbons.atoms"
        self.__rawAtomicMatrixFileName       = "rawAtomicMatrix.matrix"
        self.__filteredAtomicMatrixFileName  = "filteredAtomicMatrix.matrix"
        self.__percentageOfAlignmentFileName = "percentageFile.results"
        self.varianceFileName                = "variance.results"
        self.maskFileName                    = "varianceMask.results"
    def __load_variability_data(self, minimumPercentage):
        """
        __alphaCarbonsFileName  where to write the aligned C_alphas for each pdb.
        __alnFileName           the alignment filename. No need to be already formatted 
        __workingDirectory      where the pdbs are stored
        onlyBackbone          select only the backbone. Default False
        onlyAlphaCarbon       select only the alphaCarbons. Default True. Consider that all the latter programs work over this option and no other!!!!
        
        self.removedLines contains the removed lines, that is if a residue position is not aligned over a percentage for all
        the pdbs, this line is removed, and its order is this array        
        """
        self.__mmult = generateMultipleAlignment.getResidues(minimumPercentage, self.__alphaCarbonsFileName, self.__alnFileName, 
                                               self.__workingDirectory, False, True, False)
        
        generateMultipleAlignment.composeAtomicMatrix(self.__alphaCarbonsFileName, self.__rawAtomicMatrixFileName)
        self.removedLines = generateMultipleAlignment.filterAtomicMatrix(self.__rawAtomicMatrixFileName, self.__filteredAtomicMatrixFileName)
        self.dimensions   = generateMultipleAlignment.getDimensions(self.__filteredAtomicMatrixFileName)
        
    def __load_variability_measures(self):
        (self.__mask, self.__matrix)                    = svdAnalysis.loadCAmatrix(self.__filteredAtomicMatrixFileName, self.dimensions)
        (self.__covMatrix, self.__covMask)              = svdAnalysis.getCovarianceMatrix(self.__mask, self.__matrix)
        (self.atoms_variance, self.atoms_variance_mask) = svdAnalysis.getAtomsVariance(self.__covMatrix, self.__covMask)

    def serialize_variance_results(self):
        svdAnalysis.serializeAtomsVariance(self.atoms_variance, self.atoms_variance_mask, self.varianceFileName, self.maskFileName)
    
    def load(self, minimumPercentage, verbose):
        self.__minimumPercentage = minimumPercentage
        self.__generate_file_names()
        if verbose:
            print ("loading variability data at %.3f percentage of alignment\n" %(minimumPercentage,))
        self.__load_variability_data(minimumPercentage)
        if verbose:
            print "loading variability measures\n"
        self.__load_variability_measures()
        
        
    def loadResidues(self, maskedPDBFileName):
        """Loads for each residue in the pdb to be fitted -and contained in the multiple alignment- its corresponding position in the alignment"""
        self.__maskedPDBFileName = maskedPDBFileName         
        self.__pdbResidues = {} # key :residueInPDB, value = residueNumber

        self.__mmult.serializePercentageOfAlignment(self.__percentageOfAlignmentFileName)
        for i in range(len(self.__mmult.residues)):
            self.__pdbResidues[self.__mmult.getPositionInPDB(self.__maskedPDBFileName, i)] = i
        
    
    def compute_masked_residues(self, maskedPDBFileName, applyPercentage=False):
        """
        for each atom, give me the variance value of its residue in the multiple alignment
        stores in self.__masked_atoms a list with tuples (atom, variance,percentageOfAlignment) in order
        """
        self.__maskedPDBFileName = maskedPDBFileName
        self.__masked_atoms = []
        self.__maskedPDB = pdb.PDB(self.__maskedPDBFileName)
        self.__maskedPDB.load_atoms_hetams()
        self.__maskedPDB.loadResidues()
        self.__maskedPDB.loadIntSeqRes()
        
        atom_order = 0
        for position_in_alignment in range(len(self.__mmult.residues)):
            if self.__mmult.residues[position_in_alignment] > self.__minimumPercentage:
                positionInPDB = self.__mmult.getPositionInPDB(self.__maskedPDBFileName, position_in_alignment)
                if positionInPDB != -1:
                    atom = self.__maskedPDB.getCAlphaFromResidue(positionInPDB)
                    if atom is not None:
                        variance = self.atoms_variance[0][atom_order]
                        if applyPercentage:
                           self.__masked_atoms.append((atom, self.__mmult.residues[position_in_alignment]*variance))
                        else:
                           self.__masked_atoms.append((atom, variance))
                atom_order+=1

                    
    def serialize_masked_atoms(self, resultFileName, toBild=False):
        result = open(resultFileName, "w")
        if toBild:
            #(max, min, mean, normalizedAtoms) = svdAnalysis.atoms_variance_statistics(self.atoms_variance,self.atoms_variance_mask)
            self.computeColorsMap()
            for atom, variance in self.__masked_atoms:
                result.write(self.get_color(variance))
                result.write("\n")
                result.write(atom.toBildString())
                result.write("\n")
        else:
            for atom, variance in self.__masked_atoms:
                result.write(atom.maskedAtomToString(variance))
                result.write("\n")
                
        result.close()
    
    def generate_stats(self, bildFileName,histogramFileName):
        self.serialize_masked_atoms(resultFileName=bildFileName, toBild=True)
        self.generateHistogramOfPercentageOfAlignment(histogramFileName, show=False, serializePercentageOfAlignment=True)
    
    def computeColorsMap(self):
        differentValues = set()
        for atom, variance in self.__masked_atoms:
            differentValues.add(variance)
           
        groups = math.ceil(len(differentValues)/len(self.colors))
        l = []
        l.extend(differentValues)
        l.sort()
        i      = 0
        color  = 0
        self.__map = {}
        self.kk = self.__map
        for value in l:
            if i>groups:
               i=0
               try:
                   self.__map[value]=self.colors[color]
               except IndexError:
                    self.__map[value]=self.colors[len(self.colors)-1]
               color+=1
            else:
              i+=1
        
            
        self.__map[value]=self.colors[len(self.colors)-1]
        
        self.__maxKeys = self.__map.keys()
        self.__maxKeys.sort()

    def get_color(self, value):
        self.kkk = self.__maxKeys
        for i in range(len(self.__maxKeys)):
            if value < self.__maxKeys[i]:
                return self.__map[self.__maxKeys[i]] 
        
        return self.colors[len(self.colors)-1]
        
    def serializePercentageOfAlignment(self, percentageFileName):   
        self.__mmult.serializePercentageOfAlignment(percentageFileName)
    
    def generateHistogramOfPercentageOfAlignment(self, histogramFileName, show=False, serializePercentageOfAlignment=True):
        """note that the result file names are without extension """
        v = self.__mmult.residues.values()
        pylab.grid()
        pylab.xlabel("percentage of alignment")
        pylab.ylabel("number of aligned residues")
        #pylab.title(self.__maskedPDBFileName + (" %.3f %" % (self.__minimumPercentage,)))
        vv = pylab.hist(v, bins=len(v), normed=0)
        pylab.savefig(histogramFileName+".png")
        if show:
            pylab.show()
        if serializePercentageOfAlignment:
            self.__mmult.serializePercentageOfAlignment(histogramFileName+".per")
    

    
    