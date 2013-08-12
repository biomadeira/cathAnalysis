import time, datetime
import  MySQLdb
import re
def get_content(line):
    l = line[10:80]
    k = 10
    for i in l:
        if i == '\n':
            break
        k+=1
    return line[10:k]      
def pdbID(line):
    return get_content(line)[0:4]
def chain(line):
    return get_content(line)[4]
def start_stop(line):
        return (int(line[line.find('=')+1:line.find('STOP')]),int(line[line.find('STOP')+5:]))

def segment_delimitation_with_insertions(pattern, line):
        return re.findall(pattern, line)

def print_script(sqlScript, valueList):
    print sqlScript
    c = 0
    for i in valueList:
        print c,  i
        c+=1
    print "list's length: ", len(valueList)

def read_from_cddFile(fileName, cursor, errorFileName):
    file = open(fileName, "r")
    errorFile = open(errorFileName, "w")
    last_header = "FORMAT    "
    value = ""
    single = True

    lastDomainID = ""
    segmentIDs = []
    pattern = re.compile('([\d]+)([A-Z]*)')
    
    for line in file:
        if line.startswith("#"):
            continue
        if line.startswith("FORMAT"):
            valueList = []
            toCDDFDomains = " insert into CDDFDomains "
            toCDDFDomains+= " (FORMAT , DOMAIN , pdbID  ,chain,  VERSION, VERDATE, NAME, SOURCE , CATHCODE  , CLASS  , ARCH, TOPOL  , HOMOL  , DLENGTH, DSEQH  , DSEQS  , NSEGMENTS ) "
            toCDDFDomains+= " values (%s, %s,%s,%s,%s,%s,%s,%s,%s,%s, %s,%s,%s,%s,%s,%s,%s)"
            valueList.append(get_content(line))
            single = True
            continue
        if last_header != line[0:10] and not single:
            valueList.append(value)
            value     = ""
        last_header = line[0:10]
        
        if line.startswith("DOMAIN"):
            single = True
            valueList.append(get_content(line))
            valueList.append(pdbID(line))
            valueList.append(chain(line))
            lastDomainID = get_content(line)
            continue
        if line.startswith("VERSION"):
            single = True
            valueList.append(get_content(line))
            continue
        if line.startswith("VERDATE"):
            single = True   
            rawDate = get_content(line)
            t = time.strptime(rawDate,'%d-%b-%Y')
            dd =  datetime.date(t[0], t[1],t[2])
            valueList.append(dd)
            continue
        if line.startswith("NAME"):
            single = False   
            value+= get_content(line)
            continue
        if line.startswith("SOURCE"):    
            value+= get_content(line)
            single = False
            continue
        if line.startswith("CATHCODE"):  
            value+= get_content(line)
            single = False
            continue
        if line.startswith("CLASS"):     
            value+= get_content(line)
            single = False
            continue
        if line.startswith("ARCH"):      
            value+= get_content(line)
            single = False
            continue
        if line.startswith("TOPOL"):     
            value+= get_content(line)
            single = False
            continue
        if line.startswith("HOMOL"):     
            value+= get_content(line)
            single = False
            continue
        if line.startswith("DLENGTH"):   
            valueList.append(int(get_content(line)))
            single = True
            continue
        if line.startswith("DSEQH"):     
            value+= get_content(line)
            single = False
            continue
        if line.startswith("DSEQS"):     
            value+= get_content(line)
            single = False
            continue
        if line.startswith("NSEGMENTS"):
            single = True
            valueList.append(int(get_content(line)))
            try:
                cursor.execute(toCDDFDomains, valueList)
            except MySQLdb.Error, e:
                errorFile.write("in CDDFDomains Error %d: %s %s\n" % (e.args[0], e.args[1], lastDomainID))
            continue
        
        if line.startswith("SEGMENT"):
            valueList = []
            toSegments = "insert into CDDFSegments (segment, start, stop, startInsertion, stopInsertion, SLENGTH, SSEQH, SSEQS ) values (%s,%s,%s,%s,%s,%s,%s,%s)"
            single = True
            segmentIDs.append(get_content(line))
            valueList.append(get_content(line))
            continue
        if line.startswith("SRANGE"):
            single = True
            (start, stop, startInsertion, stopInsertion) = (-1,-1,'','')
            delimitation = segment_delimitation_with_insertions(pattern, line)
            if delimitation is not None:
                    (start, stop, startInsertion, stopInsertion) = (int(delimitation[0][0]), int(delimitation[1][0]), delimitation[0][1], delimitation[1][1])
            else:
                    errorFile.write("Error in CDDFSegments\n for line %s" % line)
            valueList.extend((start,stop,startInsertion,stopInsertion))
            continue
        if line.startswith("SLENGTH"):
            single = True
            valueList.append(get_content(line))
            continue
        if line.startswith("SSEQH"):
            single = True
            valueList.append(get_content(line))
            continue  
        if line.startswith("SSEQS"):
            value+= get_content(line)
            single = False
            continue
        if line.startswith("ENDSEG"):
            try:
                cursor.execute(toSegments, valueList)
            except MySQLdb.Error, e:
                    errorFile.write("in CDDFSegments Error %d: %s %s\n" % (e.args[0], e.args[1], lastDomainID))
            continue
        
        
        if line.startswith("//"):
            toDomainSegments = "insert into domainSegments (domain, segment) values(%s, %s)"
            for segment in segmentIDs:
                valueList = [lastDomainID, segment]
                try:
                    cursor.execute(toDomainSegments, valueList)
                except MySQLdb.Error, e:
                    errorFile.write("in domainSegments Error %d: %s %s\n" % (e.args[0], e.args[1], lastDomainID))
            lastDomainID = ""
            segmentIDs = []
            
    file.close()
    errorFile.close()
            
def operate(host, user, passwd, db, fileName, errorFileName):
    db = MySQLdb.connect(host=host,user=user, passwd=passwd,db=db)
    cursor = db.cursor()
    read_from_cddFile(fileName, cursor, errorFileName)
    cursor.close()
    db.close()
    
             
if __name__ == "__main__":
    fileName = "/home/jatienza/Desktop/cathDownloads/CathDomainDescriptionFile.v3.2.0"
    errorFileName = "/home/jatienza/Desktop/cathDownloads/CathDomainDescriptionFile.v3.2.0.errors"
    host   = 'localhost'
    user   = 'sflexfit'
    passwd = 'sflexfit'
    db     = 'sflexfit32'
    operate(host, user, passwd, db, fileName,errorFileName)