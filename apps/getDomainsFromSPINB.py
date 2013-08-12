#!/usr/bin/python
from spDomains.CATHfamilyLoader import CATHfamilyLoader
import getopt
import sys
import os
def operate(cathcode,inputDirectory="/home/inab/databases/flat/pdb", workingDirectory="/home/inab/www/CATHdomain", host="okazaki.cnb.csic.es", user="poli", passwd="R7xvgAKK", db="sflexfit", verbose=False):
    cATHfamilyLoader = CATHfamilyLoader(cathcode, host, user, passwd, db, workingDirectory)
    cATHfamilyLoader.getFamilyPDBS()
    cATHfamilyLoader.INB_downloadPDBfiles(inputDirectory, verbose=verbose)
    cATHfamilyLoader.INB_uncompressPDBFiles(inputDirectory, verbose=verbose)
    cATHfamilyLoader.extractDomains(verbose=verbose)
    compressedFileName = os.path.join(workingDirectory,cathcode)
    cATHfamilyLoader.compressResults(compressedFileName)
    logFileName = os.path.join(workingDirectory,cathcode+".log")
    #cATHfamilyLoader.deleteAll(logFileName, verbose)
    return os.path.join(workingDirectory,cathcode+".tgz")

def usage():
    print "just introduce a valid cathcode\n"

if __name__ == "__main__":
    operate(sys.argv[1])
    