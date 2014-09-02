#!/usr/bin/python

import os

def get_series_list():
    md = open("/home/mpmendenhall/Documents/aCORN/MetaData_06102014.csv").readlines()[10:]
    md = [l.split(",") for l in md]
    nprev = ""
    slist = []
    for m in md:
        if len(m)==10 and int(m[3])==1 and m[2]!=nprev:
            slist.append(int(m[2]))
            nprev = m[2]
    return slist

def check_missing():
    slist = []
    for s in get_series_list():
        fpath = os.environ["ACORN_WISHBONE"]+"/Series_%i/"%s
        if not os.path.exists(fpath):
            print "Missing",fpath
            slist.append(s)
    return slist

if __name__ == "__main__":
    
    slist = get_series_list()
    cmdlist = open("replay_cmds.txt","w")
    for s in slist:
        cmdlist.write("../../WishboneScanner %i\n"%s)
    cmdlist.close()
    
    print "Analyzing",len(slist),"series..."
    os.system("cat replay_cmds.txt")
    os.system("nice -n 15 parallel < replay_cmds.txt > replay_log.txt 2>&1")
    os.system("rm replay_cmds.txt")
    os.system("../../WishboneScanner 0");
    
    print "Checking missing files..."
    check_missing()
    print "Done."
