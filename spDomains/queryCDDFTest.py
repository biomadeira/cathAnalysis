import unittest
from spDomains.queryCDDF import *
class PDBTestCase(unittest.TestCase):
    host   = 'okazaki'
    user   = 'sflexfit'
    passwd = 'M8z2LeGR'
    db     = 'sflexfit'
    def setUp(self):
        (self.db, self.cursor) = setUpDB(PDBTestCase.host, PDBTestCase.user, PDBTestCase.passwd, PDBTestCase.db)
    def tearDown(self):
        tearDownDB(self.db, self.cursor)
        
    def test1loadCATHfamily(self):
        cathCode = '1.10.8.10'
        result = loadCATHfamily(cathCode, self.cursor)
        for pdbIDentry in result:
            print "%s" %( pdbIDentry[0],)
    
    def test2loadSegments(self):
        domainID = '1jmoL00'
        result = loadSegments(domainID, self.cursor)
        assert result is not None
        for (chain, start, stop, startInsertion, stopInsertion) in result:
            print "%s from %d to %d startInsertion[%s]  stopIsertion[%s]" % (chain, start, stop, startInsertion, stopInsertion)
            
    def test3segmentsFromPDB(self):
        pdbID =  '1jmo'
        result = segmentsFromPDB(pdbID, self.cursor)
        for (chain, start, stop, startInsertion, stopInsertion) in result:
            print "%s %d %s -- %d %s" % (chain, start,startInsertion, stop, stopInsertion) 
            
            