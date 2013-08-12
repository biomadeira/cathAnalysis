from entities.pdb import PDB
from masks.MMULT import MMULT
import libraries.OScommands
import os
def getResidues(minimumPercentage, resultsFileName, alnFileName, workingDirectory='.', onlyBackBone=False, onlyAlphaCarbon=False, writeComments = False):
    """for each pdb in the multiple structure alignment file
           for each aligned residue of the previous structure alignment file
           
                if the percentage of alignment is greater or equal than the selected percentage:
                      select every atom of the aligned residue in the structural alignment, taking into account the onlyBackBone and onlyAlphacarbon conditions
                    
     resultsFileName where to write the aligned C_alphas for each pdb.
     alnFileName  the alignment filename
     workingDirectory where the pdbs are stored
     onlyBackbone select only the backbone. Default False
     onlyAlphaCarbon select only the alphaCarbons. Default True. Consider that all the latter programs are written just considering the alpha carbons!!!!
     """
    if onlyBackBone:
        onlyAlphaCarbon = False
    if onlyAlphaCarbon:
        onlyBackBone = False
    
    workingDirectory= libraries.OScommands.gotoWorkingDirectory(True, workingDirectory)
         
    resultsFile = open(resultsFileName, "w")
    mmult = MMULT(alnFileName)
    mmult.load()
    mmult.computePercentageOfAlignment()
    mmult.loadResidues()
    
    pdbs = loadPDBS(mmult.alignments.keys(), workingDirectory)
    
    orderedPDBKeys = pdbs.keys()
    orderedPDBKeys.sort()
    
    for pdbFileName in orderedPDBKeys:
        pdb = pdbs[pdbFileName]
        if writeComments:
            resultsFile.write("#%s\n"% (pdbFileName,))
        
        orderedResiduesKeys = mmult.residues.keys()
        orderedResiduesKeys.sort()
        
        for residueNumberInAlignment in orderedResiduesKeys:
            percentage = mmult.residues[residueNumberInAlignment]
            if writeComments:
                resultsFile.write("#residueNumberInAlignment %d percentage %.3f " %(residueNumberInAlignment, percentage))
            positionInPDB = mmult.getPositionInPDB(pdbFileName, residueNumberInAlignment)
            residue = pdb.getIntSeqRes(positionInPDB)
            if residue is not None:
                if writeComments:
                    resultsFile.write(" residue in pdb %d%c \n" % (residue[0],residue[1]))
            else:
                resultsFile.write("+ + +\n")
                continue
            if percentage >= minimumPercentage:
                atoms = pdb.getAtomsFromResidue(positionInPDB)
                for atom in atoms:
                    if onlyBackBone and atom.inBackbone():
                       coordinates = atom.coordinates
                       resultsFile.write("%.3f %.3f %.3f\n" % (coordinates[0,0],coordinates[1,0], coordinates[2,0]))
                    if onlyAlphaCarbon and atom.isCAlpha():
                        coordinates = atom.coordinates
                        resultsFile.write("%.3f %.3f %.3f\n" % (coordinates[0,0],coordinates[1,0], coordinates[2,0]))
                    if not onlyBackBone and not onlyAlphaCarbon:
                       coordinates = atom.coordinates
                       resultsFile.write("%.3f %.3f %.3f\n" % (coordinates[0,0],coordinates[1,0], coordinates[2,0]))
            else:
                resultsFile.write("* * *\n") 
        resultsFile.write("#--------------------------------\n")
    
    resultsFile.close()
    return mmult

def loadPDBS(pdbFileNames, workingDirectory='./'):
    pdbs = {}
    for pdbFileName in pdbFileNames:
        pdbs[pdbFileName]= PDB(os.path.join(workingDirectory,pdbFileName))
        pdbs[pdbFileName].load_atoms_hetams()
        pdbs[pdbFileName].loadResidues()
        pdbs[pdbFileName].loadIntSeqRes()
    
    return pdbs



def composeAtomicMatrix(residuesFileName,resultsFileName, AlphaCarbon=True):
    """This function accepts a file from the former function (getResidues()). IF it is an alpha carbon there is only one line per amino acid.
      Every line of the residuesFileName can be:
      1.- atom coordinates
      2.- + + + : not aligned residue
      3.- * * * : the corresponding residue does not fit into the percentage variability
      4.- #-------: new aligned pdb
            pdb_0 pdb_1 ...pdb_j ... pdb_P
      atom_0  a_0^0          a_0^i
      atom_1
      ...
      atom_i              a_i^j
      ...
      atom_A
      
      where a_i^j = [x_i^j y_i^j z_i^j]
      """
    residuesFile = open(residuesFileName,"r")
    pdbs = []
    atoms = []
    for line in residuesFile:
        if line.startswith("+"):
            atoms.append(('+','+','+'))
            continue
        if line.startswith("*"):
            atoms.append(('*','*','*'))
            continue
        if line.startswith("#-"):#indicates the end of the atoms of the current pdb
            pdbs.append(atoms)
            atoms = []
            continue
        coordinates = line.split()
        atoms.append((float(coordinates[0]),float(coordinates[1]),float(coordinates[2])))
    residuesFile.close()
    resultsFile = open(resultsFileName,"w")
    for j in range(len(pdbs[0])):
        for i in range(len(pdbs)):
            if pdbs[i][j][0] == '+':
                resultsFile.write("%s %s %s "% ('+','+','+'))
                continue
            if pdbs[i][j][0] == '*':
                resultsFile.write("%s %s %s "% ('*','*','*'))
                continue
            resultsFile.write ("%.3f %.3f %.3f " % (pdbs[i][j][0],pdbs[i][j][1],pdbs[i][j][2]))
        resultsFile.write ("\n")
    resultsFile.close()

def filterAtomicMatrix(AtomicMatrixFileName, filteredAtomicMatrixFileName):
    """removes the lines with all * or + 
       Returns a vector with the removed positions. This is not being considered in further programs!! is a todo. Is sth that never happens ?
    """
    removedLines = []
    currentLine = 0
    AtomicMatrixFile = open(AtomicMatrixFileName, "r")
    filteredAtomicMatrixFile = open(filteredAtomicMatrixFileName, "w")
    emptyLine = True
    for line in AtomicMatrixFile:
        for coord in line.split():
            if  coord != '*' and coord != '+':
                emptyLine = False
                break
        if not emptyLine:
            filteredAtomicMatrixFile.write("%s"%(line))
        else:
            removedLines.append(currentLine)
        emptyLine = True
        currentLine+=1
        
    AtomicMatrixFile.close()
    filteredAtomicMatrixFile.close()
    return removedLines

def getDimensions(filteredAtomicMatrixFileName):
    filteredAtomicMatrixFile = open(filteredAtomicMatrixFileName, "r")
    rows = 0
    columns = 0
    columns = len(filteredAtomicMatrixFile.readline().split())
    for line in filteredAtomicMatrixFile:
        rows+=1
    rows+=1
    filteredAtomicMatrixFile.close()
    return (rows, columns)
    
def getAtomicCoordinatesMatrix(filteredAtomicMatrixFileName, finalAtomicMatrixFileName):
        """to compose the variability inter the coordinates, not to study the variability among atoms """
        filteredAtomicMatrixFile = open(filteredAtomicMatrixFileName, "r")
        finalAtomicMatrixFile    = open(finalAtomicMatrixFileName, "w")
        for line in filteredAtomicMatrixFile:
            coord = line.split()
            for x in coord[0:len(coord):3]:
                finalAtomicMatrixFile.write("%s " %(x,))
            finalAtomicMatrixFile.write("\n")
            for y in coord[1:len(coord):3]:
                finalAtomicMatrixFile.write("%s " %(y,))
            finalAtomicMatrixFile.write("\n")
            for z in coord[2:len(coord):3]:
                finalAtomicMatrixFile.write("%s " %(z,))
            finalAtomicMatrixFile.write("\n")
        filteredAtomicMatrixFile.close()
        finalAtomicMatrixFile.close()

            
                
    