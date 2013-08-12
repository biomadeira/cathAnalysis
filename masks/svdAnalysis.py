import numpy

""" """
def loadCAmatrix(caMatrixFileName, dimensions):
    """colums = number of pdbs
       rows   = 3 3D coord X CAlpha
       May include + or * that implies there is no Calpha for this atom in this pdb
     """
    caMatrixFile = open(caMatrixFileName,"r")
    
    mask     = numpy.empty(dimensions, dtype=bool)
    mask[:]  = False
    cAmatrix = numpy.empty(dimensions, dtype=float)

    row    = 0
    column = 0
    for line in caMatrixFile:
        for coord in line.split():
            if  coord != '*' and coord != '+':
                mask[row,column]=True
                cAmatrix[row,column]=float(coord)
            column+=1  
        column = 0
        row+=1                
                 
    caMatrixFile.close()
    return (mask, cAmatrix)

def serializeCAMaskAndCAMatrix(mask, cAmatrix, caMaskFileName, caMatrixFileName):
    caMaskFile   = open(caMaskFileName  , "w")
    
    for row in mask:
        for column in row:
            caMaskFile.write(str(column)) 
            caMaskFile.write(" ")
        caMaskFile.write("\n")
    
    caMaskFile.close()
    
    caMatrixFile = open(caMatrixFileName, "w")
    for row in cAmatrix:
        for column in row:
            caMatrixFile.write("%.3f"%(column,))
            caMatrixFile.write(" ")
        caMatrixFile.write("\n")
    caMatrixFile.close()

def matrixMean(mask, matrix):
    """ """
    mean       = numpy.zeros((1,3), dtype=float)
    for atomMask, atoms in zip(mask,matrix):
        validAtoms = 0
        selectedMask  = numpy.hsplit(atomMask,matrix.shape[1]/3)
        selectedAtoms = numpy.hsplit(atoms,matrix.shape[1]/3)
        for atom,atomiMask in zip(selectedAtoms, selectedMask): #for each atom
            if len(atom[atomiMask]) > 0:
                mean       += atom[atomiMask]
                validAtoms += 1
        mean/=validAtoms
        for atom,atomiMask in zip(selectedAtoms, selectedMask): #for each atom
            if len(atom[atomiMask]) > 0:
                atom-=mean[0]
    
def extractColumn(mask, matrix, column_order):
    mask_column   = mask[:,column_order]#the column_order column from mask
    matrix_column = matrix[:,column_order]#the column_order column from matrix
    return (mask_column, matrix_column) 

def extractAtomsRow(mask, matrix, i):
    return (mask[i,:], matrix[i,:])


def extractAtom(mask_column, matrix_column, row):
    atom = matrix_column[row*3:row*3+3]
    mask = mask_column  [row*3:row*3+3]
    return atom[mask]
    
def rowMeanAtom(mask,matrix):
    mean = numpy.zeros((1,3),dtype=float)
    validAtoms = 0.0
    rows = matrix.shape[0]
    selectedMask  = numpy.split(mask,rows/3)
    selectedAtoms = numpy.split(matrix,rows/3)
    for atom,atomMask in zip(selectedAtoms, selectedMask): #for each atom
        att = atom[atomMask]
        if len(att) > 0:
            mean+=att
            validAtoms+=1
    mean/=validAtoms
    return mean

def serialize_matrix(matrix, matrixFileName):
    matrixFile = open(matrixFileName, "w")
    for row in matrix:
        for f in row:
            matrixFile.write("%f\t"%(f,))
        matrixFile.write("\n")
    matrixFile.close()   
    
def getCovarianceMatrix(mask, matrix, mean=True):
    """Returns the covariance matrix and a matrix mask with the number of non valid atoms considered in the computation
       IF ONE -or more- OF THE ATOMS is not aligned, the mean of the atoms in the superfamily is considered. IF NOT:
       The couple of atoms are not considered and so the dot product sum is divided only by the number of valid atoms:
       matrix.shape[1] -1 - validAtoms
       Returns the matrix and a mask with the number of non valid atoms in each dot(row, column) operation
       Let N be the number of aligned atoms
       Let $\bar{atom_i}$ the media of the atomic coordinates of the $atom_i$ computed over all the valid aligned atoms in
       the same position in the alignment
       So, each position in the matrix is computed as follows:
       $\frac{1}{N-1}\sum_{i=1}^N \var{atom_i}\times\var{atom_i}$, with $\times$ the dot product over 3D vectors (spatial atom coordinates) 
       
    """
    matrixMean(mask, matrix)
    (rows, columns) = mask.shape
    caMatrix = numpy.zeros((rows, rows), dtype=float)
    caMask   = numpy.zeros((rows, rows), dtype=int)
    
    matrixT = matrix.T
    maskT   = mask.T
    valid = True
    dot_product = 0
    ii=0 
    for atomMask, atoms in zip(mask,matrix):#for each pdb
        selectedMask  = numpy.hsplit(atomMask,columns/3)
        selectedAtoms = numpy.hsplit(atoms,columns/3)
        jj=0
        
        for i in range(matrixT.shape[1]):
            invalidAtoms = 0
            (mask_column, matrix_column) = extractColumn(maskT, matrixT, i)
            j=0
            for atom,atomiMask in zip(selectedAtoms, selectedMask): #for each atom
                att1 = atom[atomiMask]
                att2 = extractAtom(mask_column, matrix_column, j)
                if len(att1) == 0 or len(att2) == 0:
                    valid = False
                if mean and len(att1) == 0:
                    (from_rowMask, from_row) = extractAtomsRow(mask, matrix, j) 
                    att1 = rowMeanAtom(from_rowMask, from_row)
                    invalidAtoms +=1
                    
                if mean and len(att2) == 0:
                    (from_columnMask, from_column) = extractAtomsRow(mask, matrix, i)
                    att2 = rowMeanAtom(from_columnMask, from_column)
                    invalidAtoms +=1
                if not mean and not valid:
                    pass
                else:    
                    dot_product = numpy.dot(att1,att2.T)
                j+=1
                caMatrix[ii][jj]+= dot_product
                valid       = True
                dot_product = 0
            if mean:
                caMatrix[ii][jj]/=(matrix.shape[1] -1)
            else:
                caMatrix[ii][jj]/=(matrix.shape[1] -1 - invalidAtoms)
            caMask[ii][jj] = invalidAtoms 
            jj+=1
        ii+=1
        

    return (caMatrix, caMask)


def getAtomsVariance(caMatrix, caMask):
    """
    caMatrix is a square matrix where each diagonal value is the variability associated to an aligned atom (1D variance)
    The matrix' diagonal and its corresponding mask is returned
    THIS IS NOT NECCESSARY AS THE CONVARIANCE PMATRIX  
    """
    atoms_variance      = numpy.empty((1,caMatrix.shape[0]),dtype=float)
    atoms_variance_mask = numpy.zeros((1,caMatrix.shape[0]),dtype=bool)
    for i in range(caMatrix.shape[0]):
        atoms_variance[0][i]          = caMatrix[i][i]
        if caMask[i][i] < 1:
            atoms_variance_mask[0][i] = True
            
    return (atoms_variance,atoms_variance_mask)

def serializeAtomsVariance(atoms_variance,atoms_variance_mask, varianceFileName, maskFileName):
    varianceFile = open(varianceFileName, "w")
    maskFile     = open(maskFileName, "w")
    varianceFile.write("#length %d\n"% (len(atoms_variance[0]),))
    for variance in atoms_variance[0]:
        varianceFile.write("%3.f\n"%(variance,))
    varianceFile.close()
    maskFile.write("#length %d\n"% (len(atoms_variance_mask[0]),))
    for mask in atoms_variance_mask[0]:
        maskFile.write("%d\n" % (mask,))
    maskFile.close()
    
def encarnateAtomsVariance(varianceFileName, maskFileName):
    varianceFile = open(varianceFileName,"r")
    maskFile     = open(maskFileName, "r")
    
    atoms_variance      = numpy.empty((1,int(varianceFile.readline().split()[1])),dtype=float)
    atoms_variance_mask = numpy.zeros((1,int(maskFile.readline().split()[1])),dtype=bool)
    
    i = 0
    for line in varianceFile:
        atoms_variance[0][i]=float(line)
        i+=1
    varianceFile.close()
    i = 0        
    for line in maskFile:
        atoms_variance_mask[0][i]= int(line)
    maskFile.close()
    return (atoms_variance,atoms_variance_mask)

def atoms_variance_statistics(atoms_variance,atoms_variance_mask):
    selectedAtoms = numpy.hsplit(atoms_variance, atoms_variance.shape[1])
    selectedMask  = numpy.hsplit(atoms_variance_mask, atoms_variance_mask.shape[1])
    max  = 0.0
    min  = 0.0
    mean = 0.0
    normalizedAtoms = numpy.empty((1,atoms_variance.shape[1]),dtype=float)
    for atom,atomMask in zip(selectedAtoms, selectedMask): 
        if atomMask:
            if atom > max:
                max = atom
            if atom < min:
                min = atom
            mean+=atom
    
    for i in range(atoms_variance.shape[1]):
        if atoms_variance_mask[0][i]:
            normalizedAtoms[0][i]=  atoms_variance[0][i] /mean
    
        
    return (max, min, mean, normalizedAtoms)
    