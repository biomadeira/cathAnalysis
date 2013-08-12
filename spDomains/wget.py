import os
from ftplib import FTP
def getPDBfileName(pdbCode):
    return "pdb"+pdbCode+".ent.gz"
def uncompressPDB(fileToUncompress, destFile):
    """"Uncompress a fileToUncompress file to a destFile if fileToUncompress exists. If not returns false """
    from gzip import GzipFile
    try:
        compressedPDB   = GzipFile(fileToUncompress,"r")
    except IOError:
        print "file %s does not exist" % (fileToUncompress,)
        return False
    decompressedPDB = open(destFile, "wb")
    try:
        for str in compressedPDB.read():
            decompressedPDB.write(str)
    except IOError:
        print "file %s is not a gzipped file" % (fileToUncompress,)
        return False
    
    decompressedPDB.close()
    compressedPDB.close()
    return True  

def compressFile(fileNames, compressedFileName, extension="gz"):
    """Compress a list of files
    for gzip compressed use file extension .tar.gz and modifier "w:gz"
    for bzip2 super compressed use file extension .tar.bz2 and "w:bz2"
    use "w", "w:gz" or "w:bz2" for all file types, including binary files
    http://www.daniweb.com/code/snippet630.html
    """
    
    import tarfile
    try:
        compressedPDB = tarfile.open(compressedFileName+".t"+extension, "w:"+extension)
    except IOError:
        print "file %s cannot be opened" % (compressedFileName+extension,)
        return False
    
    for fileName in fileNames:
        compressedPDB.add(fileName, os.path.basename(fileName))
    compressedPDB.close()
    return compressedFileName+".t"+extension

def setUp(site="ebi"):
    """Connects to the ftp server"""
    ftp = FTP('ftp.ebi.ac.uk')
    ftp.login("anonymous","")
    ftp.cwd('pub/databases/rcsb/pdb-remediated/data/structures/all/pdb/')
    return ftp

def tearDown(ftp):
    ftp.close()

def wgetINBPDBFileName(pdbCode,directory):
    return os.path.join(directory,"pdb"+pdbCode+".ent.Z")

def INBPDFFileExists(INBpdbFileName):
    return os.path.lexists(INBpdbFileName)

def uncompressZCATfiles(inputFileName, outputFileName):
    from popen2 import Popen4
    command = "zcat -d " + inputFileName +" > "+outputFileName
    p = Popen4(command)
    exitStatus = p.wait()
    if exitStatus != 0:
        return False
    return True


def wgetPDB(localFileName, pdbCode,ftp,site="ebi"):
    """downloads the file corresponding to the pdbCode
    If the file exists returns true, else returns false
    http://effbot.org/media/downloads/librarybook-network-protocols.pdf
    """
    value = True
    filename = getPDBfileName(pdbCode)
    file = open(localFileName, "wb")
    try:
        ftp.retrbinary('RETR ' +filename, file.write)
    except urllib2.ftplib.all_errors:
        value = False
    file.close()
    return value
