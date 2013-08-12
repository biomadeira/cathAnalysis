#!/usr/bin/python
from libraries.OScommands import executeCommand
import sys
import os
def parse_parameters(args):
    """returns an option object with the values of the optional parameters and the mandatory arguments in the list"""
    from optparse import OptionParser
    
    usage = "usage: %prog [options]  cathcode\n"
    #
    parser = OptionParser(usage=usage)
    #
    parser.add_option("-w", "--workingDirectory", dest="workingDirectory", default=".", help="directory where the resulting files will be stored") 
    parser.add_option("-m", "--mmult", dest="mmultPath", default='', help="directory of installation for mmult") 
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False, help="display program status")   
    #
    (options, args)= parser.parse_args(args)
    if len(args) != 1: #cathcode is the single mandatory parameter
        print parser.get_usage()
        sys.exit(2)
    
    return (options, args)


def operate(cathcode, workingDirectory, mmultPath, verbose):
    """TO run mmult is mandatory to move to the working directory"""
    command = os.path.join(mmultPath,"mmult ")  + cathcode + ".mmult"
    #http://docs.python.org/lib/os-process.html#os-process
    if workingDirectory == ".":
        workingDirectory = os.getcwd()
    os.chdir(workingDirectory)
    executeCommand(command, verbose)

        
    
if __name__ == "__main__":
    (optionalParameter, mandatoryParameter)=parse_parameters(sys.argv[1:])
    operate(cathcode = mandatoryParameter[0], workingDirectory = optionalParameter.workingDirectory,
             mmultPath = optionalParameter.mmultPath, verbose = optionalParameter.verbose)
    