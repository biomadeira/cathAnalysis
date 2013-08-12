from Reader import Reader
class STE (Reader):
    """"stores the steps in the execution of the whole process.
    Every state corresponds with a .py in the apps package
    All states are reset by default 
     """
     
    def __init__(self,fileName):
        self.states =  []
        Reader.__init__(self, fileName)

    
    def serialize(self,fileName):
        ste = open(fileName,"w")
        for item in self.states:
            ste.write("%s %s\n"%(item[0],str(item[1])))
        ste.close()

    def addAttributeValue(self,line):
        values = line.split()
        self.states.append((values[0], Reader.isTrue(self,values[1])))
