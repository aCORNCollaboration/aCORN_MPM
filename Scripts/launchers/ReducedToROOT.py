#!/usr/bin/python

import os

if __name__ == "__main__":
    #datadir = os.environ["ACORN_REDUCED_CALDATA"]
    datadir = os.environ["ACORN_REDUCED_DATA"]
    outdir = os.environ["ACORN_REDUCED_ROOT"]
    os.system("mkdir -p "+outdir)
    
    r2rcmd = "../../ReducedToROOT %(infl)s %(outfl)s"
    r2rZcmd = "gunzip -c %(zfl)s > %(infl)s; " + r2rcmd + "; rm %(infl)s"
    
    cmdlist = open("replay_cmds.txt","w")
    
    flist = os.listdir(datadir)
    flist.sort()
    for f in flist:
        
        if not (f[0] in ["S","s"] and ".txt" in f):
            continue
        
        if f[0]=="S": # special case for differently-capitalized filenames
            os.system("mv %s/%s %s/%s"%(datadir,f,datadir,f.lower()))
            f = f.lower()
        if "_cal" in f:
            oldf = f
            f = f.replace("_cal","")
            os.system("mv %s/%s %s/%s"%(datadir,oldf,datadir,f))
        
        outfname = "%s/%s.root"%(outdir,f.replace(".txt","").replace(".gz",""))
        if os.path.isfile(outfname):
            print outfname,"already exists. Skipping."
            continue
                
        if f[-2:] == "gz":
            zfl = datadir+"/"+f
            infl = datadir+"/tmp_"+f[:-3]
            cmdlist.write(r2rZcmd%{"zfl":zfl, "infl":infl, "outfl":outfname}+"\n")
        else:
            infl = datadir+"/"+f
            cmdlist.write(r2rcmd%{"infl":infl, "outfl":outfname}+"\n")
                
    cmdlist.close()
    os.system("cat replay_cmds.txt")
    os.system("nice -n 15 parallel < replay_cmds.txt")
    os.system("rm replay_cmds.txt")
