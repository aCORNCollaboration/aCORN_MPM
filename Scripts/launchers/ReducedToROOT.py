#!/sw/bin/python2.7

import os

if __name__ == "__main__":
	datadir = "/Users/michael/Documents/aCORN_data/"
	outdir = datadir+"/ROOTified/"
	r2rcmd = "/Users/michael/Documents/aCORN_MPM/ReducedToROOT"
	
	cmdlist = open("replay_cmds.txt","w")
		
	for d in os.listdir(datadir):
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
