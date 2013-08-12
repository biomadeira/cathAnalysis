from Reader import Reader
class CFG(Reader):
    """Reads configuration files """
    def __init__(self, fileName):
        Reader.__init__(self, fileName)
    
    def addAttributeValue(self,line):
        values = line.split(None, 1)
        self.__dict__[values[0]] = values[1]
        
