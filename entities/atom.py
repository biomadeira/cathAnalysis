from string import rjust, join
from numpy import matrix, linalg
from elements import elements
class Atom:
    """Loads an atom from a atom pdb line.
    Every operation is performed over the coordinates column matrix """
    def __init__(self, rawLine, model=0):
        if len(rawLine) < 81:
            rawLine += ' '*(81-len(rawLine))
        
        self.serial          = int(rawLine[6:11])
        self.name            = rawLine[12:16]
        self.altLoc          = rawLine[16]
        self.residue         = rawLine[17:20]
        
        self.chainIdentifier = rawLine[21]
        
        self.resSeq      = int(rawLine[22:26])
        self.iCode       =     rawLine[26] #code for insertion of residues
        
        self.coordinates = matrix([float(rawLine[30:38]),float(rawLine[38:46]), float(rawLine[46:54])] ).T

        self.occupancy = float(rawLine[54:60])
        self.tempFactor= float(rawLine[60:66])
       
        self.element         = rawLine[76:78]
        self.charge          = rawLine[78:80]
        
        self.model = model
        
    def toString(self, newSerial=0):
        sb = []
        sb.extend("ATOM  ")
        if newSerial != 0:
            sb.extend(rjust(str(newSerial),5))
        else:
            sb.extend(rjust(str(self.serial),5))
        sb.extend(' ')
        sb.extend(self.name)
        sb.extend(self.altLoc)
        sb.extend(self.residue)
        sb.extend(' ')
        sb.extend(self.chainIdentifier)
        sb.extend(rjust(str(self.resSeq),4))
        sb.extend(self.iCode)
        sb.extend('   ')
        sb.extend(rjust("%.3f" % self.coordinates[0],8))
        sb.extend(rjust("%.3f" % self.coordinates[1],8))
        sb.extend(rjust("%.3f" % self.coordinates[2],8))
        
        sb.extend(rjust("%.2f" % self.occupancy,6))
        sb.extend(rjust("%.2f" % self.tempFactor,6))
        sb.extend('          ')
        sb.extend(self.element)
        sb.extend(self.charge)
        return join(sb,'')
    
    def maskedAtomToString(self, maskValue):
        sb = []
        sb.extend("ATOM  ")
        sb.extend(rjust(str(self.serial),5))
        sb.extend(' ')
        sb.extend(self.name)
        sb.extend(self.altLoc)
        sb.extend(self.residue)
        sb.extend(' ')
        sb.extend(self.chainIdentifier)
        sb.extend(rjust(str(self.resSeq),4))
        sb.extend(self.iCode)
        sb.extend('   ')
        sb.extend(rjust("%.3f" % self.coordinates[0],8))
        sb.extend(rjust("%.3f" % self.coordinates[1],8))
        sb.extend(rjust("%.3f" % self.coordinates[2],8))
        sb.extend('   ')
        sb.extend(rjust("%.3f" % maskValue,8))
        return join(sb,'')
    
    def toBildString(self):
        sb = []
        sb.extend(".dot ")        
        sb.extend(rjust("%.3f" % self.coordinates[0],8))
        sb.extend(rjust("%.3f" % self.coordinates[1],8))
        sb.extend(rjust("%.3f" % self.coordinates[2],8))
        return join(sb,'')
    
    def coordinates(self):
        return self.coordinates
    
    def shift(self, shift):
        self.coordinates+=shift
    
    def center(self, center):
        self.coordinates-=center
    
    def rotate(self, matrix):
        """Matrix must be an square 3x3 rotatation matrix """
        self.coordinates = matrix*self.coordinates
    
    def itk_ops(self, matrix, shift):
        self.coordinates = ((self.coordinates - shift).T*matrix).T
        
    def inBackbone(self):
        if self.name.strip() == "CA" or self.name.strip() == "CB" or self.name.strip() == "O" or self.name.strip() == "N":
            return True
        return False
    def isCAlpha(self):
        if self.name.strip() == "CA":
            return True
        return False
        
    def inChain(self, chain):
        if self.chainIdentifier.strip() == chain.strip():
            return True
        return False
    def inSegment(self, start,stop):
        if self.resSeq >= start and self.resSeq <= stop:
            return True 
        return False
    
    def inInsertionsSegment(self,start, stop,startInsertion, stopInsertion, startOrder, stopOrder):
        if startInsertion == '' and stopInsertion == '':
            return self.inSegment(start, stop)
        else:           
            if self.resSeq == start:
                if startOrder:
                    if self.iCode >= startInsertion and self.resSeq >= startInsertion:
                        return True
                    else:
                        return False
                else:
                    if self.iCode <= startInsertion and self.resSeq <= startInsertion:
                        return True
                    else:
                        return False
            if self.resSeq > start and self.resSeq < stop:
                return True
            if self.resSeq == stop:
                if stopOrder:
                    if self.iCode <= stopInsertion:
                        return True
                    else:
                        return False
                else:
                    if self.iCode >= stopInsertion:
                        return True
                    else:
                        return False
            return False
    
    def atomicWeight(self):
        """Returns the atomic weight of the corresponding atom """
        return elements[self.element.strip()]
         
    def alternateLocation(self):
        return self.altLoc
    
    def insertionResidueCode(self):
        return self.resSeq
    def hasAlternateLocations(self):
        return self.altLoc != ""
    
    def inModel(self, model):
        return self.model == model 
    
    def inhighestConformation(self, anAltLoc):
        """An atom is in its highest conformation if there is no alternative location and if it's equal to anAltLoc"""
        if self.altLoc == ' ':
            return True
        else:
            if self.altLoc == anAltLoc:
                return True
        return False
    def getResidue(self):
        return (self.resSeq, self.iCode)