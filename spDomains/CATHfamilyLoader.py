import os.path
import wget
import queryCDDF 
from entities.pdb import PDB

class CATHfamilyLoader:
    """For a cath superfamily:
    1.-get every pdb that contains a domain of the superfamily from the DataBase that contains the CDDF cath file
    2.-downloads pdbs
    3.-uncompress pdbs
    4.-extracts domains
    self.pdbIDs is a tuple((pdbID,)) containing every pdb that contains a domain of the selected family
    self.notDownloadedPDBS = [] contains every pdb that does not exist in the ftp server
    self.notUncompressedFiles = [] contains every pdb that has not been uncompressed cause of it has not been downloaded or because an IOerror 
    """
    def __init__(self, cathcode, host, user, passwd, db, workingDirectory):
        self.cathcode = cathcode
        self.__host   = host
        self.__user   = user
        self.__passwd = passwd
        self.__db     = db
        self.workingDirectory = workingDirectory;
        self.notDownloadedPDBS    = []
        self.downloadedPDBS       = []
        self.notUncompressedFiles = []
        self.uncompressedFiles    = []
        self.notValidSegments = []
        self.validSegments    = [] 
    def getFamilyPDBS(self):
        """ downloads every pdb that contains a domain of the cath family  
        (a pdb can contain more than one domain of the same family)"""        
        (db, cursor) = queryCDDF.setUpDB(self.__host, self.__user, self.__passwd, self.__db)
        self.pdbIDs  = queryCDDF.loadCATHfamily(self.cathcode, cursor)
        queryCDDF.tearDownDB(db, cursor)
        
    def __composeLocalFileNameExtension(self, workingDirectory, pdbID):
        return os.path.join(workingDirectory, pdbID)
    
    def __composeLocalFileName(self, workingDirectory, fileName, extension=".pdb"):
        return os.path.join(workingDirectory, fileName+extension)
    
    def downloadPDBfiles(self, site="ebi", verbose=False):
        """ working directory must be well formed for the system the file is being executed"""
        ftp = wget.setUp()
    
        for pdbIDentry in self.pdbIDs:
            if verbose:
                print ("%s "% (pdbIDentry[0],))
            if not wget.wgetPDB(self.__composeLocalFileNameExtension(self.workingDirectory, wget.getPDBfileName(pdbIDentry[0])), pdbIDentry[0], ftp, site):
                self.notDownloadedPDBS.append(wget.getPDBfileName(pdbIDentry[0]))
                if verbose:
                    print (" not downloaded \n")
            else:
                self.downloadedPDBS.append(wget.getPDBfileName(pdbIDentry[0]))
                if verbose:
                    print (" downloaded \n")
        
        wget.tearDown(ftp)
        
    def INB_downloadPDBfiles(self,directory="/home/inab/databases/flat/pdb", verbose=False):
        for pdbIDentry in self.pdbIDs:
            if verbose:
                print ("%s "% (pdbIDentry[0],))
            if not wget.INBPDFFileExists(wget.wgetINBPDBFileName(pdbIDentry[0], directory)):
                    self.notDownloadedPDBS.append(wget.wgetINBPDBFileName(pdbIDentry[0], directory))
                    if verbose:
                        print (" not downloaded \n")
            else:
                    self.downloadedPDBS.append(wget.wgetINBPDBFileName(pdbIDentry[0], directory))
                    if verbose:
                        print (" downloaded \n")
        
    def uncompressPDBFiles(self, verbose=False):
        for pdbIDentry in self.pdbIDs:
            if verbose:
                print ("%s "% (pdbIDentry[0],))
            if not wget.uncompressPDB(self.__composeLocalFileNameExtension(self.workingDirectory, wget.getPDBfileName(pdbIDentry[0])), \
                          self.__composeLocalFileName(self.workingDirectory, pdbIDentry[0])):
                self.notUncompressedFiles.append(self.__composeLocalFileName(self.workingDirectory, pdbIDentry[0]))
                if verbose:
                    print (" not uncompressed \n")
            else:
                self.uncompressedFiles.append(self.__composeLocalFileName(self.workingDirectory, pdbIDentry[0]))
                if verbose:
                    print (" uncompressed\n")
    
    def INB_uncompressPDBFiles(self, inputDirectory="/home/inab/databases/flat/pdb", verbose=False):
        for pdbIDentry in self.pdbIDs:
            if verbose:
                print ("%s "% (pdbIDentry[0],))
            if not wget.uncompressZCATfiles(wget.wgetINBPDBFileName(pdbIDentry[0],inputDirectory), os.path.join(self.workingDirectory,pdbIDentry[0]+".pdb")):
                self.notUncompressedFiles.append(os.path.join(self.workingDirectory,pdbIDentry[0]+".pdb"))
                if verbose:
                    print (" not uncompressed \n")
            else:
                self.uncompressedFiles.append(os.path.join(self.workingDirectory,pdbIDentry[0]+".pdb"))
                if verbose:
                    print (" uncompressed\n")
    
    def extractDomains(self, model=0, onlyBackbone=False, verbose=False):
        """1.-query the database to know every pair(pdbID, domainID)
           2.-for each domain of the family contained in the pdb, load its segments.
           3.-for each domain serialize it with the name domain.pdb
        """
        (db, cursor)      = queryCDDF.setUpDB(self.__host, self.__user, self.__passwd, self.__db)
        self.PDB_Domains  = queryCDDF.loadCATHfamilyDomains(self.cathcode, cursor)
        for (pdbID, domain) in self.PDB_Domains:
            validDomain =  True
            segments    = queryCDDF.loadSegments(domain, cursor)
            for segment in segments:
                if segment[1] != -1:
                    self.validSegments.append((pdbID, domain, segment))
                else:
                    self.notValidSegments.append((pdbID, domain, segment))
                    validDomain = False
                    
            if validDomain:
                try:
                    pdb = PDB(self.__composeLocalFileName(self.workingDirectory, pdbID, ".pdb"))
                    pdb.load_atoms_hetams()
                except IOError:
                    print "%s does not exist" %(self.__composeLocalFileName(self.workingDirectory, pdbID, ".pdb"),)
                else:
                    pdb.serializeSegments(segments, self.__composeLocalFileName(self.workingDirectory, domain, ".pdb"), model, onlyBackbone)
                    if verbose:
                        print ("%s domain extracted" % (domain,))
                    
        queryCDDF.tearDownDB(db, cursor)
        
    def deleteCompressedPDBS(self,verbose=False):
        """Delete compressed pdb files, raw pdbs. Just keep the extracted domains and the print logFile"""
        for pdbID in self.downloadedPDBS:
            fileName = self.__composeLocalFileNameExtension(self.workingDirectory, pdbID)
            try:
                os.remove(fileName)
                if verbose:
                    print("removing %s\n"% (fileName,))
            except OSError:
                print("%s not found \n"% (fileName,))
        
        for pdbID in self.pdbIDs:
            fileName = self.__composeLocalFileName(self.workingDirectory, pdbID[0], ".pdb")
            try:
                os.remove(fileName)
                if verbose:
                    print("removing %s\n"% (fileName,))
            except OSError:
                print("%s not found \n"% (fileName,))
        
        self.printResults()
                
    def deleteAll(self, logFileName=None,verbose=False):
        """delete the pdbs, the extracted domains and the log file """
        for pdbID in self.downloadedPDBS:
            fileName = self.__composeLocalFileNameExtension(self.workingDirectory, pdbID)
            try:
                os.remove(fileName)
                if verbose:
                    print("removing %s\n"% (fileName,))
            except OSError:
                print("%s not found \n"% (fileName,))
        #
        for pdbID in self.pdbIDs:
            fileName = self.__composeLocalFileName(self.workingDirectory, pdbID[0], ".pdb")
            try:
                os.remove(fileName)
                if verbose:
                    print("removing %s\n"% (fileName,))
            except OSError:
                print("%s not found \n"% (fileName,))
        #       
        domains = set()
        for (pdbID, domain) in self.PDB_Domains:
            domains.add(domain)
        for domain in domains:
            domainFileName = self.__composeLocalFileName(self.workingDirectory, domain, ".pdb") 
            os.remove(domainFileName)
            if verbose:
                print ("%s deleted\n"%(domainFileName,))
        #
        if not logFileName is None:
            os.remove(logFileName)
            if verbose:
                print ("%s deleted\n"% (logFileName,))
        else:
            os.remove(os.path.join(self.workingDirectory,self.cathcode+".log"))
            if verbose:
                print ("%s deleted\n"% (os.path.join(self.workingDirectory,self.cathcode+".log"),))
        
        os.remove(os.path.join(self.workingDirectory,self.cathcode+".mmult"))
        if verbose:
                print ("%s deleted\n"% (os.path.join(self.workingDirectory,self.cathcode+".mmult"),))
    
    
    def compressResults(self, compressedFileName, extension="gz"):
        """Compress the domain files and the logFile. Returns the tgz file """
        fileNames = []
        domains = set()
        for (pdbID, domain) in self.PDB_Domains:
            domains.add(domain)
        for domain in domains:
            fileNames.append(self.__composeLocalFileName(self.workingDirectory, domain, ".pdb"))
        fileNames.append(self.printResults())
        fileNames.append(os.path.join(self.workingDirectory,self.cathcode+".mmult"))
        return wget.compressFile(fileNames, compressedFileName, extension)
        
     
    def generateDomainList(self,verbose=False):
        file = open(os.path.join(self.workingDirectory,self.cathcode+".mmult"),"w")
        file.write("MAMMOTH\n")
        domains = set()
        for (pdbID, domain) in self.PDB_Domains:
            domains.add(domain)
        for domain in domains:
            file.write(domain+".pdb\n")
        file.close()
        if verbose:
                print ("%s domain list generated \n"%(os.path.join(self.workingDirectory,self.cathcode+".mmult"),))            
        
    def printResults(self, fileName=None):
        """Somehow a fileName must be substituted. I'm duplicating code! """
        if fileName is None:
            fileName = os.path.join(self.workingDirectory,self.cathcode+".log")
        file = open(fileName, "w")
        file.write("NOT downloaded pdbs \n")
        for pdbID in self.notDownloadedPDBS:
            file.write("%s \n" %(pdbID,))
        
        file.write("downloaded pdbs \n")
        for pdbID in self.downloadedPDBS:
            file.write("%s\n" %(pdbID,))
        
        file.write("NOT uncompressed files \n")
        for fileID in self.notUncompressedFiles:
            file.write("%s\n" %(fileID,))
                
        file.write("uncompressed files \n")
        for fileID in self.uncompressedFiles:
            file.write("%s\n" %(fileID,))
                
        file.write("NOT valid segments\n")
        for pdbID, domain, segment in self.notValidSegments:
            file.write("pdb[%s]_domain[%s]_chain[%s]_start[%d] stop[%d]\n" % (pdbID, domain, segment[0], segment[1], segment[2]))
            
        file.write("valid segments \n")
        for pdbID, domain, segment in self.validSegments:
            file.write("pdb[%s]_domain[%s]_chain[%s]_start[%d] stop[%d]\n" % (pdbID, domain, segment[0], segment[1], segment[2]))
            
        file.write("generated domains \n")
        for pdbID, domain, segment in set(self.validSegments):
            file.write("%s.pdb\n" % (domain,))
            
        return fileName
                
        
        