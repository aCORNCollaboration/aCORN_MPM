#!/usr/bin/python

import os

if __name__ == "__main__":
    datadir = os.environ["ACORN_REDUCED_DATA"]
    outdir = os.environ["ACORN_REDUCED_ROOT"]
    os.system("mkdir -p "+outdir)
    r2rcmd = "../../ReducedToROOT"
    
    cmdlist = open("replay_cmds.txt","w")
    
    flist = os.listdir(datadir)
    flist.sort()
    for f in flist:
        if f[0] in ["S","s"] and f[-4:] == ".txt":
            if f[0]=="S": # special case for differently-capitalized filenames
                os.system("mv %s/%s %s/%s"%(datadir,f,datadir,f.lower()))
                f = f.lower()
            #if f[1:5] not in ["1519"]:
            #    continue
            cmdlist.write(r2rcmd + " %s/%s"%(datadir,f) + " %s/%s.root\n"%(outdir,f[:-4]))
    cmdlist.close()
    os.system("cat replay_cmds.txt")
    os.system("nice -n 15 parallel < replay_cmds.txt")
    os.system("rm replay_cmds.txt")
