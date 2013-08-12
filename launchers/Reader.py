class Reader:
    def __init__(self, fileName):
        file = open(fileName, "r")
        i = 0
        for line in file:
            if self.__isFormattable(line):
                try:
                    self.addAttributeValue(line.strip())
                except BaseException , Exception :
                    print ("Error formatting line [%d] \n \t%s\n"%(i,line))
            i+=1
        file.close()
        
    def __isFormattable(self,line):
        tLine = line.strip()
        if tLine.startswith("#"):
            return False
        if len(tLine) < 1:
            return False
        return True
        
    
    def isTrue(self, string):
        string = string.lower()
        if string == "true" or string == 't':
            return True
        return False
    
    def format(self, type, value):
        """is called in the constructor """
        if type == 'str':
            return value
        if type == 'bool':
            return self.isTrue(value)
        if type == 'int':
            return int(value)
        if type == 'float':
            return float(value)
        if type == 'scp':
            return (value,'scp')
        if type == 'cmd':
            return (value,'cmd')
    
    def __getParameterName(self,parameter):
        return parameter.lstrip("-")
    
    def command(self,line):
        values = line.split()
        script = [values[0],' ']
        extending = False
        
        for value in values[1:]:
            if value == '+':
                extending = True
                continue
            if hasattr(self,value):
                if not extending:
                    script.append(' ')
                script.append(self.__dict__[value])
                extending = False
                continue
            if value.startswith("\""):
                value=value.lstrip('"').rstrip('"')
                if not extending:
                    script.append(' ')
                script.append(value)
                extending = False
                continue
            if value.startswith("--"):
                script.append(' ')
                self.appendParameter(script, value)
                extending = False


        return (''.join(script),'cmd') #returns an iterable to set the script names apart
    
    def appendParameter(self, script, value):
        #if the parameter needs a parameter value
        if type(True) != type(self.__dict__[self.__getParameterName(value)]):
            script.append(value)
            script.append(' ')
            script.append(str(self.__dict__[self.__getParameterName(value)]))
        else:##if the parameters are boolean, that its, they don't need a parameter value ie --verbose 
            if self.__dict__[self.__getParameterName(value)]:
                script.append(value)
                script.append(' ')
                #script.append(str(self.__dict__[self.__getParameterName(value)]))
                    
    def appendMandatoryArgument(self,script,value):
        if hasattr(self, value):
            script.append(self.__dict__[value])
        else:
            script.append(value)
                    
    def script(self, line):
        """set up the parameters of the scripts """
        values = line.split()
        script = []
        for value in values:
            if value.startswith("--"): 
                self.appendParameter(script, value)
            else:#is a mandatory argument (without --)
                self.appendMandatoryArgument(script, value)
            script.append(' ')
    
        return (''.join(script),'scp') #returns an iterable to set the script names apart