from string import rjust, join
from numpy import matrix, linalg
from atom import Atom
class Hetatm(Atom):
    def init(self, rawString,modelNumber):
        Atom.init(self, rawString,modelNumber)
    def toString(self):
        sb = []
        sb.extend("HETATM")
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