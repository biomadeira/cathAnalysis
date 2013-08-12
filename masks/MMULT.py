class MMULT:
    """This class loads a multiple alignment file (clustal like) from mmult (mammoth multiple alignment)
    to use this class:
    1.-load the alignment file (self.alignments is populated)
    2.-compute the percentage of alignment 
    Afterwards, you can extract the residue number of an aligned pdb from the position it occupies 
    """
    def __init__(self, alnFileName):
        self.__file = open(alnFileName,"r")
        self.alignments = {} # key:pdb file name value: set of aminoacids and ---
        self.lenOfAlignment = -1
        self.residues = {}
        self.pdbResiduesOrder = {}
    
    
    def load(self):
        for line in self.__file:
            if not line.isspace() and not line.startswith("CLUSTAL") and not line.startswith(" "):
                fLine = line.split()
                if self.alignments.has_key(fLine[0]):
                    previous = self.alignments[fLine[0]]
                    self.alignments[fLine[0]] = previous +fLine[1]
                else:
                   self.alignments[fLine[0]] = fLine[1]
        
        self.lenOfAlignment = len(self.alignments[fLine[0]])
        self.__file.close()
        
    def serialize(self, fileName):
        file = open(fileName, "w")
        keys = self.alignments.keys()
        keys.sort()
        for key in keys:
            file.write("%s %s \n" % (key, self.alignments[key]))
        file.close()
        
    def computePercentageOfAlignment(self):
        for i in range(self.lenOfAlignment):
            self.residues[i]= 0.0
        for alignment in self.alignments.items():
            i = 0
            for residue in alignment[1]:
                if residue != '-':
                    self.residues[i]+=1
                i+=1
        
        for i in self.residues.keys():
            value = float(self.residues[i])
            self.residues[i] = value /len(self.alignments)
    
    def findFirstAlignedResidue(self, pdbFileName):
        """Returns the first aligned residue in the structural alignment regarding the first aligned residue of the pdbs considered"""
        i = 0
        for residue in self.alignments[pdbFileName]:
            if residue == '-':
                i+=1
            else:
                return i
    
    def loadResidues(self):
        for pdbFileName in self.alignments.keys():
            self.pdbResiduesOrder[pdbFileName]=self.findFirstAlignedResidue(pdbFileName)

    def getPositionInPDB(self, pdbFileName, residueNumberInAlignment):
        """returns the position occupied in the pdb -residueNumber- of the residue in the alignment -residueNumberInAlignment-
           If residueNumberInAlignment > len(self.alignments[pdbFileName]) it returns the position of the last position in 
           the alignment"""
        alignedResidues    = 0
        notAlignedResidues = 0
        for residue in self.alignments[pdbFileName]:
            if residue != '-':
                alignedResidues+=1
            else:
                notAlignedResidues+=1
            if (alignedResidues + notAlignedResidues) > residueNumberInAlignment:
                break
        if residue != '-':
            return alignedResidues - 1
        else:        
            return -1
        
    def serializePercentageOfAlignment(self, percentageFileName):
        """Serialize the percentage of alignment of the aligned residues"""
        percentageFile = open(percentageFileName,"w")
        orderedKeys = self.residues.keys()
        orderedKeys.sort()
        for residueNumberInAlignment in orderedKeys:
            percentageFile.write(" %d %.3f \n" % (residueNumberInAlignment, self.residues[residueNumberInAlignment]))
        percentageFile.close()
        
