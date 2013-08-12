def executeCommand(command, verbose=True):
    from popen2 import Popen4
    if verbose:
        from os import getcwd
        print command
        print getcwd()
    p = Popen4(command)
    exitStatus = p.wait()
    if exitStatus != 0:
        raise OSError, command +" fail error: " + str(exitStatus) 



def gotoWorkingDirectory(workingDirectoryExists=False, workingDirectory='.'):
    import os
    workingDirectory=workingDirectory.strip()
    if workingDirectoryExists:
        if workingDirectory == ".":
            workingDirectory = os.getcwd()
        if workingDirectory[0]=="/":
            pass
        if workingDirectory[0:2]=="./":
            workingDirectory = os.path.join(os.getcwd(),workingDirectory[2:])
        if workingDirectory[0:2]=="..": #relative path ../
            os.chdir(workingDirectory) 
            workingDirectory = os.getcwd()
    else:
        workingDirectory = os.getcwd()
    os.chdir(workingDirectory)
    return workingDirectory

def splitFileName(fileName):
    """capable of managing fileNames with multiple points n.a.m.e.extension -> (n.a.m.e, extension) """
    b = os.path.basename(fileName).split(".")
    sb = []
    for i in range(len(b)-1):
        sb.append(b[i])
        sb.append('.')
    try:
        sb.pop()
    except IndexError:
        return (b,"")
    return ("".join(sb), b[len(b)-1])


def isTrue(string):
    f = string.strip().lower()()
    if f == "true" or f == 't':
        return True
    return False 