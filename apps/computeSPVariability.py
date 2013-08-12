#!/usr/bin/python
from masks.GenerateMask import MaskDataLoader

def operate(workingDirectory, alnFileName, minimumPercentage,stats, patternPDBFileName, bildFileName, histogramFileName, verbose):
    dataLoader = MaskDataLoader(workingDirectory, alnFileName, patternPDBFileName, bildFileName, histogramFileName)
    dataLoader.operate(minimumPercentage, stats, verbose)


def parse_parameters(args):
    """returns an option object with the values of the optional parameters and the mandatory arguments in the list"""
    #http://docs.python.org/library/optparse.html
    from optparse import OptionParser, OptionGroup
    usage = "Compute the SP variability from a multiple alignment file in CLUSTAL format.\n" 
    usage+= "usage: %prog alignment file [options] \n"
    usage+= "Example:\n"
    usage+= "./computeSPVariability 1.10.100.180\n"
    
    usage+= "the generated files are:\n"
    usage+= "alphaCarbons.atoms\n"
    usage+= "rawAtomicMatrix.matrix\n"
    usage+= "filteredAtomicMatrix.matrix\n"
    usage+= "percentageFile.results\n"
    usage+= "variance.results\n"
    usage+= "varianceMask.results\n"
    
    parser = OptionParser(usage=usage)
    parser.add_option("-p","--percentage",dest="percentage",default=0.2,help="minimum percentage of alignment of the selected residues")
    parser.add_option("-v","--verbose",action="store_true", dest="verbose",default=True,help="display program status")   
    
    parser.add_option("-w","--workingDirectory", dest="workingDirectory",default=".",help="directory where the  files will be stored") 
    
    group = OptionGroup(parser, "Statistical Information", "statistical file names")
    group.add_option("-s","--stats",action="store_true", dest="stats",default=False,help="Generate statistical information about the superfamily")
    group.add_option("-t","--pattern", dest="patternPDBFileName", default="",help="PDB file to perform the analisis")       
    group.add_option("-b","--bild", dest="bildFileName", default="bild.bild",help="File containing the atoms grouped by variability")
    group.add_option("-g","--histogram", dest="histogramFileName", default="histogram",help="File containing the histrogram file")
    parser.add_option_group(group)
    us
    (options, args)= parser.parse_args(args)

    if len(args) != 1: #cathcode is the single mandatory parameter
        print parser.get_usage()
        sys.exit(2)
    
    return (options,args)

if __name__ == "__main__":
    import sys
    (optionalParameter,mandatoryParameter)=parse_parameters(sys.argv[1:])    
    operate(workingDirectory=optionalParameter.workingDirectory, alnFileName=mandatoryParameter[0],
                 minimumPercentage=float(optionalParameter.percentage),stats= optionalParameter.stats,
                 patternPDBFileName=optionalParameter.patternPDBFileName, 
                 bildFileName= optionalParameter.bildFileName, histogramFileName=optionalParameter.histogramFileName,
                verbose = optionalParameter.verbose)

    