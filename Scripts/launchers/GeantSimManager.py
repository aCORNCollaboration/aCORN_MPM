#!/usr/bin/python
import os
import time
from math import *
from optparse import OptionParser

# killall -9 GeantSimManager.py; killall -9 parallel; killall -9 aCORN_G4_prod

class GeantSimManager:
	
	def __init__(self, simName, vacuum="1.e-7 torr", fmap=None):
		
		self.settings = {}
		
		self.settings["run_num"] = 0
		self.settings["simName"] = simName
		self.settings["vacuum"] = vacuum
		self.settings["gunpos_mm"] = [0.,0.,0.]
		self.settings["fieldmapcmd"] = "#/detector/fieldmapfile UNUSED"
		if fmap:
			self.settings["fieldmapcmd"] = "/detector/fieldmapfile "+fmap
		self.settings["physlist"] = "livermore"
		self.settings["ana_args"] = ""
		self.settings["extra_cmds"] = ""
		self.settings["gen_cmds"] = ""
		self.settings["vis_cmd"] = ""
		
		self.g4_out_dir_base = None
		
		self.anagroup = 10 # number of files to group together for final analyzer result
				
	def enable_vis(self):
		# open selected visualization driver
		self.settings["vis_cmd"] = "/vis/open OGLIX\n"

		# set viewpoint
		self.settings["vis_cmd"] += "/vis/viewer/set/viewpointThetaPhi 105 0\n"
		self.settings["vis_cmd"] += "/vis/viewer/panTo -0.1 0\n"
		#self.settings["vis_cmd"] += "/vis/viewer/zoom 0.75\n"
		
		# enable drawing detector geometry
		self.settings["vis_cmd"] += "/vis/drawVolume\n"
		
		# draw more complex wireframe
		#self.settings["vis_cmd"] += "/vis/viewer/set/auxiliaryEdge true\n"
		# solid surface drawing options
		#self.settings["vis_cmd"] += "/vis/viewer/set/style surface\n"
		#self.settings["vis_cmd"] += "/vis/viewer/set/hiddenEdge 1"
		
		self.settings["vis_cmd"] += "/tracking/verbose 2\n"
		
		# track visualization
		if 1:
			self.settings["vis_cmd"] += "/vis/modeling/trajectories/create/drawByCharge myTrackVis\n"
			# extra tracking points
			self.settings["vis_cmd"] += "/vis/modeling/trajectories/myTrackVis/default/setDrawStepPts true\n"
			self.settings["vis_cmd"] += "/vis/modeling/trajectories/myTrackVis/default/setDrawAuxPts true\n"
			
			self.settings["vis_cmd"] += "/vis/modeling/trajectories/select myTrackVis\n"
			self.settings["vis_cmd"] += "/vis/scene/add/trajectories\n"
			self.settings["vis_cmd"] += "/vis/scene/add/trajectories rich\n"
			self.settings["vis_cmd"] += "/vis/scene/add/hits\n"
		
		self.settings["vis_cmd"] += "/vis/viewer/flush\n"
	
	def set_evtsrc(self,evtsrc):
		self.settings["evtsrc"] = evtsrc;
		
	def set_dirs(self):
		self.g4_workdir = os.environ["G4WORKDIR"]
		self.g4_bindir = os.environ["G4BINDIR"]
		self.type_dir = self.settings["simName"]
		if "evtsrc" in self.settings:
			self.g4_evtsdir = os.environ["G4EVTDIR"]+"/"+self.settings["evtsrc"]
			self.type_dir +="_"+self.settings["evtsrc"]
		
		if "gunenergy" in self.settings:
			self.type_dir += "_%.1fkeV"%self.settings["gunenergy"]
			self.settings["gen_cmds"] += "/gun/energy %f keV"%self.settings["gunenergy"]
		
		if not self.g4_out_dir_base:
			self.g4_out_dir_base = self.g4_workdir+"/output/"
		self.g4_out_dir = self.g4_out_dir_base+"/%s/"%self.type_dir
		self.g4_log_dir = self.g4_workdir+"/logs/%s/"%self.type_dir
		self.g4_macro_dir = self.g4_workdir+"/macros/%s/"%self.type_dir
		self.g4_out_name = "%s/g4_run_%%s.root"%self.g4_out_dir

	def launch_sims(self,maxIn=100000,hours_old=0):
		
		self.set_dirs()
		parallel_jobfile = "%s/jobs.txt"%self.g4_macro_dir
		
		os.system("mkdir -p %s"%self.g4_macro_dir)
		os.system("mkdir -p %s"%self.g4_out_dir)
		os.system("mkdir -p %s"%self.g4_log_dir)
		
		nruns = 1
		if "evtsrc" in self.settings:
			inflist = [f for f in os.listdir(self.g4_evtsdir) if f[:5]=="Evts_"]
			inflist.sort()
			inflist = inflist[:maxIn]
			nruns = len(inflist)
		
		oldtime = time.time() - hours_old*3600
		
		# main simulations
		os.system("rm -r %s/*"%self.g4_macro_dir)
		jobsout = open(parallel_jobfile,"w")
		ucnG4_prod = self.g4_bindir+"/aCORN_G4_prod"
		onejob = ""
		
		# set up macros for each job
		for rn in range(nruns):
			self.settings["run_num"] += 1
			self.settings["jobname"] = self.settings["simName"]+"_%i"%self.settings["run_num"]
			self.settings["outfile"]=self.g4_out_name%str(self.settings["run_num"])
			if "evtsrc" in self.settings:
				self.settings["evtfile"]="/benchmark/gun/evtfile "+self.g4_evtsdir+"/"+inflist[rn]
			else:
				self.settings["evtfile"]=""
			self.settings["nevt"] = min(maxIn,10000) # assume this many events per input file... TODO something more elegant
			self.settings["joblog"] = "%s/gen_macro_%i.txt"%(self.g4_log_dir,self.settings["run_num"])
			g4_sub_file = "%s/geantjob_%i.sub"%(self.g4_macro_dir,self.settings["run_num"])
			
			self.settings["gunpos"] = "%g %g %g mm"%tuple(self.settings["gunpos_mm"])
			
			# skip recently-run jobs
			if os.path.exists(self.g4_out_name%str(self.settings["run_num"])) and os.stat(self.g4_out_name%str(self.settings["run_num"])).st_mtime > oldtime:
				continue;
			
			# generate macro file
			open(os.path.expanduser("%s/geantgen_%i.mac"%(self.g4_macro_dir,self.settings["run_num"])),"w").write(open("GeantGenMacroTemplate.mac","r").read()%self.settings)
			# single job execution command, appended to batch job file
			onejob = ucnG4_prod + " %s/geantgen_%i.mac %s"%(self.g4_macro_dir,self.settings["run_num"],self.settings["physlist"])
			jobsout.write(onejob+" > %s 2>&1\n"%self.settings["joblog"])
		
		jobsout.close()
		
		print "Running simulation jobs..."
		os.system("cat "+parallel_jobfile)
		if nruns > 1:
			os.system("nice -n 20 parallel < %s"%parallel_jobfile)
		else:
			os.system(onejob)
		os.system("rm "+parallel_jobfile)
	
	
	def launch_postanalyzer(self,nMin=0,nMax=100000):
		print "Running post analyzer..."
		self.set_dirs()
		resim_jobfile = "%s/resim_jobs.txt"%self.g4_macro_dir
		jobsout = open(resim_jobfile,"w")
		anafiles = [ (int(f[:-5].split("_")[-1]),self.g4_out_dir+"/"+f) for f in os.listdir(self.g4_out_dir) if f[:7]=="g4_run_"]
		anafiles.sort()
		nanalyzed = 0
		self.settings["analyzer"]="UCNA_MC_Analyzer"
		while anafiles:
			outlist_name = self.g4_out_dir+"outlist_%i.txt"%nanalyzed
			fout = open(outlist_name,"w")
			for f in anafiles[:self.anagroup]:
				fout.write(f[1]+"\n")
			fout.close()
			anafiles = anafiles[self.anagroup:]
			print "\n----- %s ------"%outlist_name
			os.system("cat "+outlist_name)
			analyzer_bin = self.g4_bindir+"/"+self.settings["analyzer"]
			if nMin <= nanalyzed <= nMax:
				jobsout.write("%s %s %s/analyzed_%i.root%s\n"%(analyzer_bin,outlist_name,self.g4_out_dir,nanalyzed,self.settings["ana_args"]))
			nanalyzed += 1
		jobsout.close()
		print "\n----- %s ------"%resim_jobfile
		os.system("cat "+resim_jobfile)
		print
		os.system("nice -n 10 parallel < %s"%resim_jobfile)
		os.system("rm %s/outlist_*.txt"%self.g4_out_dir)
		os.system("rm "+resim_jobfile)


if __name__ == "__main__":
	
	parser = OptionParser()
	parser.add_option("-k", "--kill", dest="kill", action="store_true", default=False, help="kill running replays")
	options, args = parser.parse_args()
	if options.kill:
		os.system("killall -9 parallel")
		os.system("killall -9 ucnG4_prod")
		os.system("killall -9 UCNA_MC_Analyzer")
		os.system("killall -9 GeantSimManager.py")
		exit(0)
	
	# self.settings["ana_args"] += " saveall"

	####################
	# beta decays
	####################

	# unpolarized beta baseline: 5e7 in 520 clusters
	if 0:
		betaSim = GeantSimManager("20120823")
		betaSim.settings["physlist"]="livermore"
		betaSim.set_evtsrc("neutronBetaUnpol")
		betaSim.g4_out_dir_base = "/data2/mmendenhall/G4Out/2010/"
		betaSim.launch_sims(nEvents=5e7,nClusters=520,hours_old=10*24)
		betaSim.launch_postanalyzer()
		exit(0)


	####################
	# calibration sources
	####################


	####################				
	# scintillator electron gun
	####################
	if 0:
		for l in [100, 200, 400, 600, 800]:
			iline = GeantSimManager("eGun")
			iline.set_generator("eGunRandMomentum")
			iline.settings["positioner"] = "Fixed"
			iline.settings["gunenergy"] = l
			iline.settings["ana_args"] += " undead saveall"
			iline.launch_sims(nEvents=1e6,nClusters=36,hours_old=0)
			iline.launch_postanalyzer()
		

	##################
	# visualization test
	##################
	if 1:
		vtest = GeantSimManager("Vis_Test")
		vtest.settings["gunenergy"] = 500
		vtest.enable_vis()
		vtest.launch_sims(maxIn=100)


