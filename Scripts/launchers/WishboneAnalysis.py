#!/usr/bin/python

import os

if __name__ == "__main__":

    md = open("/home/mpmendenhall/Documents/aCORN/MetaData_06102014.csv").readlines()[10:]
    md = [l.split(",") for l in md]
    
    cmdlist = open("replay_cmds.txt","w")
    nprev = ""
    for m in md:
        if len(m)==10 and int(m[3])==1 and m[2]!=nprev:
            cmdlist.write("../../WishboneScanner %i\n"%int(m[2]))
            nprev = m[2]
    cmdlist.close()
    
    os.system("cat replay_cmds.txt")
    os.system("nice -n 15 parallel < replay_cmds.txt > replay_log.txt 2>&1")
    os.system("rm replay_cmds.txt")
