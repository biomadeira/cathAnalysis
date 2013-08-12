from entities.atom import Atom
from entities.hetatm import Hetatm
from entities.aminoacids import aminoacids

from numpy import matrix, dot, sqrt
from libraries.rotationMatrices import *
class PDB:
    """Please, take care of working with different models, I haven't checked it properly
    For each entity in the PDB, the model number is stored
     """
    def __init__(self, pdbFileName):
        self.__file = open(pdbFileName, "r")
        self.__fileName = pdbFileName
        self.atoms = {}
        self.hetatms = {}
        self.__rest = {}
        self.__numberOfLines = 0
        self.__firstResidue = 0
        
        self.residues  = {}
        self.intSeqRes = {}
        
        self.CAlphas = {}

    def load_atoms_hetams(self):
        """Loads every atom and hetam to its corresponding structure. The key of each line is its line number 
        Each atom/hetam holds the model it belongs to.
        rest is for when only rotation or translation operations are performed.
        """
        i = 0
        modelNumber = 0
        for line in self.__file:
            if line.startswith("ENDMDL"):
                modelNumber+=1
                self.__rest[i] = line
            else:
                if line.startswith("ATOM"):
                   self.atoms[i] = Atom(line,modelNumber)
                else:
                   if line.startswith("HETATM"):
                        self.hetatms[i]= Hetatm(line,modelNumber)
                   else:
                        if not line.isspace():
                             self.__rest[i] = line
            i+=1
        self.__numberOfLines = i
        self.__file.close()
        
        #To extract the atoms of a residue
        if len(self.atoms)>0:
            firstAtomOrder = min(self.atoms.keys())
            self.__firstResidue = (self.atoms[firstAtomOrder].resSeq,self.atoms[firstAtomOrder].iCode) 
        #print "atoms " + str(len(self.atoms)) 
        #print "hetatms " + str(len(self.hetatms))
        #print "rest " + str(len(self.__rest))
        
    def serialize(self, newPDBfileName):
        """This method doesn't consider the model """
        newFile = open(newPDBfileName, "w")
        i = 0
        while i < self.__numberOfLines:
            if self.atoms.get(i):
                newFile.write(self.atoms.get(i).toString())
                newFile.write("\n")
            else:
                if self.hetatms.get(i):
                    newFile.write(self.hetatms.get(i).toString())
                    newFile.write("\n")
                else:
                    if self.__rest.get(i):
                        newFile.write(self.__rest.get(i))
            i+=1
            
        newFile.close()
    
    def shift(self, vector):
        """Adds vector to every atom and hetatm coordinates"""
        for a in self.atoms.values():
            a.shift(vector)
        for h in self.hetatms.values():
            h.shift(vector)
    
    def rotate(self, rotationMatrix, center):
        """ Rotates every atom considering center as the center of rotation"""
        for a in self.atoms.values():
            a.center(center)
            a.rotate(rotationMatrix)
        for h in self.hetatms.values():
            h.center(center)
            h.rotate(rotationMatrix)
    
    def chimera_transform(self, center, rotationMatrix, shift):
        """For each atom, hetatm is first centered in the center of coordinates, then rotated, then shifted, and finally recentered
        transform = lambda atom, center, shift, rotationMatrix: rotationMatrix*(atom - center)+shift+center
        """
        for a in self.atoms.values():
            a.center(center)
            a.rotate(rotationMatrix)
            a.shift(shift+center)
        for h in self.hetatms.values():
            h.center(center)
            h.rotate(rotationMatrix)
            h.shift(shift+center)
    def ash_transform(self, center, rotationMatrix, shift):
        """For each atom, hetatm is first centered in the center of coordinates, then rotated an finally shifted
        transform = lambda atom, center, shift, rotationMatrix: rotationMatrix*(atom - center)+shift
        """
        for a in self.atoms.values():
            a.center(center)
            a.rotate(rotationMatrix)
            a.shift(shift)
        for h in self.hetatms.values():
            h.center(center)
            h.rotate(rotationMatrix)
            h.shift(shift)
            
    def itk_transform(self, matrix, shift):
            for a in self.atoms.values():
                a.itk_ops(matrix, shift)
            for h in self.hetatms.values():
                h.itk_ops(matrix, shift)
    
    def getHighestOccupancyLocation(self, anAtom):
        """Some of the atoms of the residue can have an unique conformation (== ' ') and so an occupancy factor or 1
        This algorithm returns the location with highest occupancy 
         """
        occupancy = anAtom.occupancy
        location  = anAtom.altLoc
        for atom in self.atoms.values():
            if anAtom.model == atom.model:
                if anAtom.resSeq == atom.resSeq:
                    if atom.altLoc != anAtom.altLoc:
                        if atom.altLoc != ' ' and atom.occupancy > anAtom.occupancy:
                            occupancy = atom.occupancy
                            location = atom.altLoc
        return location
            
    def getAtomsFromSegment(self,chain, start, stop, startInsertion, stopInsertion, startOrder, stopOrder, model=0, onlyBackbone=False):
        """ returns a dictionary of the atoms of the backbone in a chain in a segment
        Please, take into account that the atoms are NOT stored in ORDER. Those are in a dictionary 
        The atoms that not fulfil the conditions are not included in the dictionary
        """
        segmentAtoms = {}
        highestAltLocation = ''
        
        for (lineNumber,atom) in self.atoms.items():
            if atom.inModel(model):
                if atom.inChain(chain):
                    if atom.inInsertionsSegment(start, stop,startInsertion, stopInsertion, startOrder, stopOrder):
                        if onlyBackbone and atom.inBackbone():
                              if atom.hasAlternateLocations():
                                  highestAltLocation = self.getHighestOccupancyLocation(atom)
                                  if atom.inhighestConformation(highestAltLocation):
                                      segmentAtoms[lineNumber] = atom
                              else:
                                  segmentAtoms[lineNumber] = atom
                        else: 
                              if atom.hasAlternateLocations():
                                  highestAltLocation = self.getHighestOccupancyLocation(atom)
                                  if atom.inhighestConformation(highestAltLocation):
                                      segmentAtoms[lineNumber] = atom
                              else:
                                  segmentAtoms[lineNumber] = atom
                            
        return segmentAtoms
    
    def getAtomsFromDomain(self, segments, model=0, onlyBackbone=False):
        """returns a set of the atoms of the backbone or all in a chain for every segment of the domain
        segments = ((chain, start, stop, startInsertion, stopInsertion))
         """
        domainAtoms = []
        for segment in segments:
            startOrder = self.orderOfInsertions(segment[1])
            stopOrder  = self.orderOfInsertions(segment[2])
            domainAtoms.append(self.getAtomsFromSegment(segment[0], int(segment[1]),int(segment[2]),segment[3],segment[4],startOrder, stopOrder, model, onlyBackbone))
        return domainAtoms
    
    def getRestForSegments(self, segments):
        """@TODO Get the rest of the structures  to  compose a pdb """
        pass
    
    def serializeSegments(self, segments, domainFileName, model=0,onlyBackbone=False):
        """Stores the selected segments of the pdb in the domainFileName.
        The segments structure is a list -segments = []- composed by -segments.append() - a (chain, start, stop) segments
        A new serial for each atom is set. If more than 1 model is found, the first model is used. 
         """
        file = open(domainFileName, "w")
        domains = self.getAtomsFromDomain(segments,model,onlyBackbone)
        newSerial = 1
        for segment in domains:
            orderedKeys = segment.keys()
            orderedKeys.sort()
            for lineNumber in orderedKeys:
                file.write("%s\n" %(segment[lineNumber].toString(newSerial),))
                newSerial+=1
        file.write("%s\n" % ("END"+77*" ", ))
        file.close()
        
    def centerOfMass(self,model=0):
        """This method considers only the atoms in a model"""
        cm     = matrix([0.0, 0.0, 0.0]).T
        weight = 0.0
        for a in self.atoms.values():
            if a.inModel(model):
                cm    += a.atomicWeight()*a.coordinates
                weight+= a.atomicWeight()
        return cm/weight
    
    
    def sphereRadius(self,model=0):
        """returns the radius of the sphere the pdb is contained, and the center of mass """
        cm = self.centerOfMass(model) 
        radius = 0.0
        for a in self.atoms.values():
            if a.inModel(model):
                dist_vector = (a.coordinates - cm).A.ravel()
                distance    = sqrt(dot(dist_vector,dist_vector))
                print distance
                if distance > radius:
                    radius = distance
        return (cm, radius)
    
         
    
    def orderOfInsertions(self, resSeq):
        """The order of the insertions is defined as the backwards or forwards alphasbetical order that is A,B,C ... or ...C, B, A
        to do this, a couple of consecutive atoms sharing the order, with an iCode != '' are found, and after that, the order is determined.
        True if the order is from A to Z False otherwise.
        """
        firstAtom = None
        for atom in self.atoms.values():
            if atom.resSeq == resSeq and  atom.iCode != '' and firstAtom is  None:
                firstAtom = atom
                continue
            if firstAtom is not None and atom.resSeq == resSeq and atom.iCode !='' and firstAtom.iCode != atom.iCode:
                return firstAtom.iCode < atom.iCode
            else:
                firstAtom = None
                
        return True
    def fasta(self, fastaFileName, model=0):
        """generates a fasta file of the pdb selecting the aminoacids of the atoms in the model and /or in the backbone """
        fastaFile = open(fastaFileName,"w")
        fastaFile.write(">%s Model %d  \n" % (self.__fileName, model))
        keys = self.atoms.keys()
        keys.sort()
        resSeq  = -1
        iCode = ''
        currentLine = []
        for line in keys: 
            if self.atoms[line].inModel(0):
                if self.atoms[line].resSeq != resSeq or self.atoms[line].iCode != iCode:
                    if len(currentLine) < 79:
                       currentLine.append(aminoacids[self.atoms[line].residue])
                    else:
                        currentLine.append(aminoacids[self.atoms[line].residue])  
                        fastaFile.write("%s\n" % ''.join(currentLine))
                        currentLine = []
                    resSeq = self.atoms[line].resSeq
                    iCode = self.atoms[line].iCode
        fastaFile.write("%s\n" % ''.join(currentLine))
        
        fastaFile.close()
        
    def getAtomsFromResidue(self, intSeqRes):
        """"Returns the atoms associated with a residue. The intResSeq parameter is the order the residue occupies in its pdb,
        considering the iCodes (insertion codes).LoadResidues and LoadIntSeqRes  methods must be called before using getAtomsFromResidue"""
        atoms = []
        residue = self.getIntSeqRes(intSeqRes)
        for atom in self.atoms.values():
            if atom.getResidue() == residue:
                atoms.append(atom)
        return atoms
    
    def getCAlphaFromResidue(self, intSeqRes):
        """Returns the CAlpha of the residue.The intResSeq parameter is the order the residue occupies in its pdb,
        considering the iCodes (insertion codes).LoadResidues and LoadIntSeqRes  methods must be called before using getAtomsFromResidue"""
        residue = self.getIntSeqRes(intSeqRes)
        for atom in self.atoms.values():
            if atom.isCAlpha() and atom.getResidue() == residue:
                return atom
    
    def getFirstResidue(self):
        """returns the first residue in the pdb defined as the (atom.SeqRes,atom.iCode)  """
        return self.__firstResidue
    
    def getIntSeqRes(self, intSeqRes):
        """Returns the residue defined as (atom.seqRes, atom.iCode) """
        return self.intSeqRes.get(intSeqRes)
    
    def loadIntSeqRes(self):
        """to work,  loadResidues must be called before"""
        for (key, value) in self.residues.items():
            self.intSeqRes[value] = key
    
    def loadResidues(self):
        """loads the residues of the pdb, defining a residue as (atom.seqRes, atom.iCode) for any different combination in the pdb """
        intResSeq = 0
        for atom in self.atoms.values():
            if  not self.residues.has_key(atom.getResidue()):
                self.residues[atom.getResidue()]= intResSeq
                intResSeq+=1
                
    
        