#!python

import os

if __name__ == "__main__":
    datadir = os.env["ACORN_REDUCED_DATA"]
    outdir = os.env["ACORN_REDUCED_ROOT"]
    os.system("mkdir -p "+outdir)
    r2rcmd = "../../ReducedToROOT"
    
    cmdlist = open("replay_cmds.txt","w")
    
    for d in os.listdir(datadir):
        #if d[:3] not in ["cal"]:
        #continue
        subdirs = []
        try:
            subdirs = os.listdir(datadir+"/"+d)
            except:
                continue
            for f in subdirs:
                if f[-4:] == ".txt":
                    cmdlist.write(r2rcmd + " %s/%s/%s"%(datadir,d,f) + " %s/%s.root\n"%(outdir,f[:-4]))
                    cmdlist.close()
                    os.system("cat replay_cmds.txt")
                    os.system("nice -n 15 parallel < replay_cmds.txt")
                    os.system("rm replay_cmds.txt")
