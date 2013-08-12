#!/usr/bin/python
import sys
import os
from Reader import Reader
from cfgReader import CFG
from steReader import STE
from libraries.OScommands import executeCommand
class Launcher(Reader):
    
    """Contains the name, type and default value of the parameters """
    def __init__(self, lchFileName):
        self.__scripts = dict()
        self.__cmds    = dict()
        Reader.__init__(self, lchFileName)
            
    
    def addAttributeValue(self,line):
        """Instantiates the values for the program and the states
        attName   type  value (delimited by blank spaces)"""
        values = line.split(None, 2)
        att = Reader.format(self, values[1], values[2])
        if hasattr(att, '__iter__'):
            if att[1] == 'scp':
                self.__scripts[values[0]] = att[0]
            if att[1] == 'cmd':
                self.__cmds[values[0]] = att[0]
        else:
            self.__dict__[values[0]] = att
        
    def operate(self,verbose=True):
        """
        1.-loads CFG and STE files
        2.-uses CFG file
        3.-executes scripts
        4.-serialize states
         """
        self.__cfg = CFG(self.cfg)
        self.__ste   = STE(self.ste)


        """ import os
        if hasattr(self, 'workingDirectory'):
            if workingDirectory == ".":
                self.workingDirectory = os.getcwd()
            os.chdir(workingDirectory)
        else:
            os.chdir(os.getcwd())"""
        
        #load configuration parameters. 
        for (name, value) in self.__cfg.__dict__.items():
            if type(self.__dict__[name]) == type(True):
                self.__dict__[name]=Reader.isTrue(self,value)
                continue
            if type(self.__dict__[name]) == type(float):
                self.__dict__[name]=float(value)
                continue
            if type(self.__dict__[name]) == type(int):
                self.__dict__[name]=int(value)
                continue
            else:
                self.__dict__[name]=value
        
        #compose the scripts
        for (name, script) in self.__scripts.items():
            self.__scripts[name] = Reader.script(self, script)[0]
        for (name, cmd) in self.__cmds.items():
            self.__cmds[name] = Reader.command(self, cmd)[0]
        
        #execute the scripts corresponding to the enabled states

        i = 0    
        for (name, state) in self.__ste.states:
            if not state :
                print ("%s\n"% (name,))
                if self.__scripts.has_key(name):
                    executeCommand(self.__scripts[name])
                else:
                    if self.__cmds.has_key(name):
                        import os
                        """1.-move to the working directory """
                        if self.__dict__.has_key('workingDirectory'):
                            if self.__dict__['workingDirectory'] == ".":
                                os.chdir(os.getcwd())
                            else:
                                os.chdir(self.__dict__['workingDirectory'])
                        else:
                            os.chdir(os.getcwd())
                        executeCommand(self.__cmds[name])
                    else:
                        raise KeyError("script not command found "+name)
                        
                                
                self.__ste.states[i]=(name,True)
            i+=1    
        self.__ste.serialize(self.ste)

    

def usage(self):
    print("Usage: \n")
    print("./launcher.py <fileName>.lch\n")

def filter(args):
    if len(args) > 1:
        usage()
        sys.exit()

if __name__ == "__main__":
    filter(sys.argv[1:])
    launcher = Launcher(sys.argv[1])
    launcher.operate()