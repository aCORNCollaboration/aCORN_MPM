#include "aCORN_MC_Analyzer.hh"
#include <TMath.h>
#include <cfloat>
#include <string.h>

aCORN_MC_Analyzer::aCORN_MC_Analyzer(const std::string& outfname): aCORN_G4_Analyzer(outfname), saveAllEvents(false) { }

void aCORN_MC_Analyzer::setupOutputTree() {
	printf("Adding branches for aCORN_MC_Analyzer...\n");
	anaTree->Branch("Edep",&Edep,"Edep/D");
	anaTree->Branch("EdepQ",&EdepQ,"EdepQ/D");
	anaTree->Branch("EdepAll",&EdepAll,"EdepAll/D");
	
	char tmp[1024];
	sprintf(tmp,"EdepSD[%i]/D",N_SD);
	anaTree->Branch("EdepSD",EdepSD,tmp);
	sprintf(tmp,"thetaInSD[%i]/D",N_SD);
	anaTree->Branch("thetaInSD",thetaInSD,tmp);
	sprintf(tmp,"thetaOutSD[%i]/D",N_SD);
	anaTree->Branch("thetaOutSD",thetaOutSD,tmp);
	sprintf(tmp,"keInSD[%i]/D",N_SD);
	anaTree->Branch("keInSD",keInSD,tmp);
	sprintf(tmp,"keOutSD[%i]/D",N_SD);
	anaTree->Branch("keOutSD",keOutSD,tmp);
	sprintf(tmp,"hitCountSD[%i]/I",N_SD);
	anaTree->Branch("hitCountSD",hitCountSD,tmp);
	
	anaTree->Branch("ScintPos",ScintPos,"ScintPos[3]/D");
	anaTree->Branch("ScintPosSigma",ScintPosSigma,"ScintPosSigma[3]/D");
}

void aCORN_MC_Analyzer::resetAnaEvt() {
	Edep = EdepQ = 0;
	for(AxisDirection d=X_DIRECTION; d<=Z_DIRECTION; ++d)
		ScintPos[d] = ScintPosSigma[d] = 0;
	EdepAll = 0;
	for(size_t ii=0; ii<N_SD; ii++) {
		hitCountSD[ii] = keInSD[ii] = keOutSD[ii] = thetaInSD[ii] = thetaOutSD[ii] = EdepSD[ii] = 0.;
		hitTimeSD[ii]=FLT_MAX;
	}
}

void aCORN_MC_Analyzer::processTrack() {
	
	// detector ID numbers
	const int ID_scint = 0;
	
	// particle ID numbers
	const int PDG_ELECTRON = 11;
	const int PDG_POSITRON = -11;
	
	const double kMe = 510.998910; // electron mass, keV
	
	// total deposited energy in all sensitive volumes
	EdepAll += trackinfo->Edep;
	if(detectorID < N_SD) {
		EdepSD[detectorID] += trackinfo->Edep;
		if(trackinfo->isEntering && pID==PDG_ELECTRON)
			hitCountSD[detectorID]++;
		if(pID==PDG_ELECTRON && trackinfo->hitTime<hitTimeSD[detectorID]) {
			hitTimeSD[detectorID] = trackinfo->hitTime;
			double pin_x = trackinfo->pIn[0];
			double pin_y = trackinfo->pIn[1];
			double pin_z = trackinfo->pIn[2];
			double pout_x = trackinfo->pOut[0];
			double pout_y = trackinfo->pOut[1];
			double pout_z = trackinfo->pOut[2];
			if(pin_x || pin_y || pin_z) {
				double magpin2 = pin_x*pin_x+pin_y*pin_y+pin_z*pin_z;
				keInSD[detectorID] = sqrt(magpin2+kMe*kMe)-kMe;
				thetaInSD[detectorID] = acos(fabs(pin_z/sqrt(magpin2)));
			}
			if(pout_x || pout_y || pout_z) {
				double magpout2 = pout_x*pout_x+pout_y*pout_y+pout_z*pout_z;
				keOutSD[detectorID] = sqrt(magpout2+kMe*kMe)-kMe;
				thetaOutSD[detectorID] = acos(fabs(pout_z/sqrt(magpout2)));
			}
		}
	}
	
	// scintillator deposited energy, position, hit time
	if(detectorID==ID_scint[s]) {
		Edep += trackinfo->Edep;
		EdepQ += trackinfo->EdepQuenched;
		for(AxisDirection d=X_DIRECTION; d<=Z_DIRECTION; ++d) {
			ScintPos[d] += trackinfo->edepPos[d];
			ScintPosSigma[d] += trackinfo->edepPos2[d];
		}
	}
}

void aCORN_MC_Analyzer::processEvent() {
	// normalize position variables
	for(AxisDirection d=X_DIRECTION; d<=Z_DIRECTION; ++d) {
		if(Edep > 0) {
			ScintPos[d] /= Edep;
			ScintPosSigma[d] = sqrt(ScintPosSigma[d]/Edep-ScintPos[d]*ScintPos[d]);
		}
	}
}

int main(int argc, char** argv) {
	
	if(argc<3) {
		cout<<"Syntax: "<<argv[0]<<" <filename containing list of raw root files> <output root file name> [saveall|undead]"<<endl;
		exit(1);
	}
	
	aCORN_MC_Analyzer UMA(argv[2]);
	
	for(int i=3; i<argc; i++) {
		std::string arg = argv[i];
		if(arg=="saveall") {
			UMA.saveAllEvents = true;
		} else {
			cout<<"Unknown argument: "<<argv[0]<<endl;
			exit(1);
		}
	}
	
	UMA.analyzeFileList(argv[1]);
	
	return 0;
}
