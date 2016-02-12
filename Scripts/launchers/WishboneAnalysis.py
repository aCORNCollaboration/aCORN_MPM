#!/usr/bin/python

import os
from optparse import OptionParser

def get_series_list_NG6():
    """NG-6 series list from metadata"""
    md = open("/home/mpmendenhall/Documents/aCORN/Reference/MetaData_06102014.csv").readlines()[10:]
    md = [l.split(",") for l in md]
    nprev = ""
    slist = []
    for m in md:
        if len(m)==10 and int(m[3])==1 and m[2]!=nprev:
            slist.append(int(m[2]))
            nprev = m[2]
    slist.sort()
    return slist

def get_series_list_analyzed():
    """Series list based on analyzed files"""
    slist = [int(f.split("_")[-1]) for f in os.listdir(os.environ["ACORN_WISHBONE"]) if f[:7] == "Series_"]
    slist.sort()
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
    
    parser = OptionParser()
    parser.add_option("--replot", dest="replot", action="store_true", default=False, help="Update analysis plots")
    parser.add_option("--analyze", dest="analyze", action="store_true", default=False, help="Analyze all data")    
    options, args = parser.parse_args()
    
    cmdbase = "../../WishboneScanner -%i\n"
    njobs = 4
    if options.analyze:
        options.replot = False
        njobs = 1
        cmdbase = "../../WishboneScanner %i\n"
    elif not options.replot:
        exit(0)
            
    slist = get_series_list_analyzed()
    cmdlist = open("replay_cmds.txt","w")
    for s in slist:
        cmdlist.write(cmdbase%s)
    cmdlist.close()
    
    print "Analyzing",len(slist),"series..."
    os.system("cat replay_cmds.txt")
    os.system("nice -n 15 parallel -j %i < replay_cmds.txt > replay_log.txt 2>&1"%njobs)
    os.system("rm replay_cmds.txt")
    os.system("../../WishboneScanner 0");
    
    print "Checking missing files..."
    check_missing()
    print "Done."
