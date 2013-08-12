#!/usr/bin/python
from spDomains.CATHfamilyLoader import CATHfamilyLoader
import sys
from os import path

def operate(cathcode, workingDirectory, host, user, passwd, db, verbose,tar,justTar):
    if verbose:
        print ("getting domains for %s superfamily"% (cathcode,))
    cATHfamilyLoader = CATHfamilyLoader(cathcode, host, user, passwd, db, workingDirectory)
    if verbose:
        print "getting information from DB"
    cATHfamilyLoader.getFamilyPDBS()
    if verbose:
        print "downloading pdb files"
    cATHfamilyLoader.downloadPDBfiles(verbose=verbose)
    if verbose:
        print "uncompressing pdb files"
    cATHfamilyLoader.uncompressPDBFiles(verbose=verbose)
    if verbose:
        print "extracting domains"
    cATHfamilyLoader.extractDomains(verbose=verbose)
    if verbose:
        print "generating domain list for mmult"
    cATHfamilyLoader.generateDomainList(verbose=verbose)
    #
    compressedFileName = path.join(workingDirectory,cathcode)
    #
    if tar:
        cATHfamilyLoader.compressResults(compressedFileName)
        logFileName = path.join(workingDirectory,cathcode+".log")
        if justTar:
             cATHfamilyLoader.deleteAll(logFileName, verbose=True)
        if verbose:
           print ("finished. You can find the resulting tgz at %s" % (path.join(workingDirectory,cathcode+".tgz"),))

    else:
        cATHfamilyLoader.deleteCompressedPDBS(verbose)
    
       
    if verbose and not tar:
        print ("finished. You can find the log file at %s" % (path.join(workingDirectory,cathcode+".log"),))
    

def parse_parameters(args):
    """returns an option object with the values of the optional parameters and the mandatory arguments in the list"""
    from optparse import OptionParser, OptionGroup
    
    usage = "usage: %prog [options]  cathcode\n"
    usage+= "Example: get data for the 1.10.10.180 cath superfamily\n"
    usage+= "./getDomainsFromSP.py  1.10.10.180\n"
    usage+= "Example:get data for the 1.10.10.180 cath superfamily with all posible conf parameters \n"
    usage+= "./getDomainsFromSP.py  1.10.10.180 -w . -m localhost -u sflexfit -p sflexfit -d sflexfit\n"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-t","--tar",dest="tar",action="store_true",default=False,help="generate a compressed tar file with the results")
    parser.add_option("-w","--workingDirectory", dest="workingDirectory",default=".",help="directory where the resulting files will be stored") 
    parser.add_option("-v","--verbose",action="store_true", dest="verbose",default=True,help="display program status")   
    parser.add_option("-j","--justTar",dest="justTar",action="store_true",default=False,help="remove all the generated files leaving just the tar")   
    
    group = OptionGroup(parser, "DataBase options", "DataBase connection parameters")
    group.add_option("-m","--DBServer", dest="host", default="okazaki.cnb.csic.es",help="database server")
    group.add_option("-u","--user", dest="user", default="poli",help="database user name")
    group.add_option("-p","--passwd", dest="passwd", default="R7xvgAKK",help="user password")
    group.add_option("-d","--db", dest="db", default="sflexfit",help="database name")
    parser.add_option_group(group)
    
    (options, args)= parser.parse_args(args)

    if len(args) != 1: #cathcode is the single mandatory parameter
        print parser.get_usage()
        sys.exit(2)
    
    return (options,args)

    
if __name__ == "__main__":
    (optionalParameter,mandatoryParameter)=parse_parameters(sys.argv[1:])
    
    #cathcode, workingDirectory, host, user, passwd, db, verbose,tar
    operate(cathcode=mandatoryParameter[0],workingDirectory=optionalParameter.workingDirectory,host=optionalParameter.host,
            user=optionalParameter.user,passwd=optionalParameter.passwd,db=optionalParameter.db,
            verbose=optionalParameter.verbose,tar=optionalParameter.tar,justTar=optionalParameter.justTar)
    