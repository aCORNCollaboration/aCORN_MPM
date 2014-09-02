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
        if f[0]=="s" and f[-7:] == "rd2.txt": # and f[-4:] == ".txt":
            cmdlist.write(r2rcmd + " %s/%s"%(datadir,f) + " %s/%s.root\n"%(outdir,f[:-4]))
    cmdlist.close()
    os.system("cat replay_cmds.txt")
    os.system("nice -n 15 parallel < replay_cmds.txt")
    os.system("rm replay_cmds.txt")
