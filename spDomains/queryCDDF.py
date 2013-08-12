import  MySQLdb

def setUpDB(host, user, passwd, db):
    db = MySQLdb.connect(host=host,user=user, passwd=passwd,db=db)
    cursor = db.cursor()
    return (db, cursor)

def tearDownDB(db, cursor):
    cursor.close()
    db.close()    

def execute_query(cursor, sqlQuery):
    cursor.execute(sqlQuery)
    return cursor.fetchall()

def execute_script(cursor, sqlScript, valuesList):
    cursor.execute(sqlScript,valuesList)

def loadCATHfamily(cathCode, cursor):
    """returns every domain as a part of every pdb in the cath superfamily"""
    query = " select distinct pdbID"
    query+= " from CDDFDomains "
    query+= " where cathcode = '"
    query+= cathCode
    query+= "'"
    query+= " order by pdbID"
    cursor.execute(query)    
    return cursor.fetchall()

def loadCATHfamilyDomains(cathCode, cursor):
    """returns every domain as a part of every pdb in the cath superfamily"""
    query = " select pdbID, DOMAIN"
    query+= " from CDDFDomains "
    query+= " where cathcode = '"
    query+= cathCode
    query+= "'"
    query+= " order by pdbID, DOMAIN"
    cursor.execute(query)    
    return cursor.fetchall()

def loadSegments(domainID, cursor):
    """returns start and stop amino acids for every segment of the domain"""
    query = " select d.chain, segs.start, segs.stop, segs.startInsertion, segs.stopInsertion"
    query+= " from CDDFDomains d, domainSegments ds, CDDFSegments segs "
    query+= "where d.domain = '"
    query+= domainID
    query+= " ' "
    query+= " and d.DOMAIN = ds.domain "
    query+= " and ds.segment = segs.segment "
    cursor.execute(query)
    return cursor.fetchall()

def domainsFromPDB(pdbID, cursor):
    """returns domains of a certain pdb"""
    query = " select domain "
    query+= " from CDDFDomains " 
    query+= " where pdbID = '"
    query+= pdbID 
    query+= " '"
    cursor.execute(query)
    return cursor.fetchall()




def segmentsFromPDB(pdbID, cursor):
    domains = domainsFromPDB(pdbID, cursor)
    segments = []
    for domain in domains:
         result = loadSegments(domain[0], cursor)
         for (chain, start, stop, startInsertion, stopInsertion) in result:
             segments.append((chain, start, stop,startInsertion, stopInsertion))
    return segments
    