#!/usr/bin/python
from masks.GenerateMask import MaskGenerator
import getopt
import sys

def operate(workingDirectory, alnFileName, minimumPercentage, pdbFileName, resultPDBFileName,verbose):
    varianceFileName = "variance.results"
    maskFileName     = "varianceMask.results"
    maskGenerator    = MaskGenerator(workingDirectory, alnFileName, minimumPercentage, 
                                  varianceFileName, maskFileName, pdbFileName, resultPDBFileName)
    maskGenerator.operate(verbose)
    
def parse_parameters(args):
    """returns an option object with the values of the optional parameters and the mandatory arguments in the list"""
    from optparse import OptionParser
    
    usage = "This program generates a mask over a pdb in the superfamily whose variability data have\n"
    usage+= "been previously computed. The order of the mandatory parameters: alignment file & pdbFile is mandatory\n"
    usage+= "usage: %prog  alignmentFile pdbFile [options] \n"
    usage+= "example:\n"
    usage+= "./generateMask alignmentFile.aln 1AH02.pdb -w .  -p 0.2  -r 1AH02_mask.pdb"
    
    parser = OptionParser(usage=usage)
    
    parser.add_option("-w","--workingDirectory", dest="workingDirectory",default=".",help="directory where the resulting files will be ")
    parser.add_option("-a","--alignmentFile", dest="alnFileName",default=".",help="alignment file")
    parser.add_option("-p","--percentage", dest="minimumPercentage",default=0.2,help="percentage of alignment. Must be equal to the percentage used to compute the SP variability")
    parser.add_option("-r","--maskResult", dest="resultPDBFileName",default="results.pdb",help="file name of the pdbresult")
    parser.add_option("-v","--verbose",action="store_true", dest="verbose",default=True,help="display program status")   
    (options, args)= parser.parse_args(args)

    if len(args) != 2: #cathcode is the single mandatory parameter
        print parser.get_usage()
        sys.exit(2)
    
    return (options,args)

if __name__ == "__main__":
    (optionalParameter,mandatoryParameter)=parse_parameters(sys.argv[1:])
    operate(workingDirectory=optionalParameter.workingDirectory, alnFileName=mandatoryParameter[0], minimumPercentage=float(optionalParameter.minimumPercentage), pdbFileName=mandatoryParameter[1], resultPDBFileName=optionalParameter.resultPDBFileName,verbose=optionalParameter.verbose)
