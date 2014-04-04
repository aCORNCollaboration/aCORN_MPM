#include "Field.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>
#include <fstream>
#include <string>
#include <cassert>

Field::Field(const G4String& filename): addAFP(false), rmax2(20*20*cm2), fieldScale(1.0) {
	LoadFieldMap(filename);
}

void Field::LoadFieldMap(const G4String& filename) {
	
	Bpoints.clear();
	Zpoints.clear();
	
	if(filename==""){
		//default field profile
		addPoint(-3.0*m,0.6*tesla);
		addPoint(-2.2*m,0.6*tesla);
		addPoint(-1.5*m,1.0*tesla);
		addPoint(1.5*m,1.0*tesla);
		addPoint(2.2*m,0.6*tesla);
		addPoint(3.0*m,0.6*tesla);
	} else {
		// load profile from file
		ifstream fin;
		fin.open(filename.data());
		if(!fin) {
			G4cout << "Can not open " << filename << G4endl;
			exit(1);
		} 
		G4String stmp;
		assert(false); // TODO new string classes
		//stmp.ReadLine(fin);  //skip the first title line
		
		std::string stmp1, stmp2;
		while(fin){
			fin >> stmp1 >> stmp2;
			//if(stmp1.IsNull()) break;
			assert(false); // TODO new string classes
			//else {
				G4double z = atof(stmp1.c_str())*m;
				G4double B = atof(stmp2.c_str())*tesla;
				addPoint(z,B);
				cout << z/m << " " << B/tesla << endl;
			//}
		}
		fin.close();
	}
}

void Field::GetFieldValue(const G4double Point[3], G4double *Bfield) const {
	
	G4double z=Point[2];	// point z
	unsigned int zindex = int(lower_bound(Zpoints.begin(), Zpoints.end(), z)-Zpoints.begin());	// location in points list
	
	// no field defined outside experimental volume
	if(zindex==0 || zindex>=Zpoints.size() || Point[0]*Point[0]+Point[1]*Point[1]>rmax2 || !fieldScale) {
		Bfield[0] = Bfield[1] = Bfield[2] = 0;
		return;
	}
	
	// interpolate between defined regions
	G4double base = 0.5*(Bpoints[zindex-1]+Bpoints[zindex]);// midpoint value
	G4double amp = 0.5*(Bpoints[zindex-1]-Bpoints[zindex]);	// variation amplitude between ends
	G4double dz = Zpoints[zindex]-Zpoints[zindex-1];		// z distance between ends
	G4double l = (z-Zpoints[zindex-1])/dz;					// fractional distance between ends
	
	Bfield[2] = base*fieldScale;
	if(amp) {
		Bfield[2] += amp*cos(l*M_PI)*fieldScale; // interpolate B_z component with cosine
		// B_r component to obey Maxwell equation grad dot B = dB_z/dz + 1/r d(r B_r)/dr = 0
		G4double Brtemp = amp*M_PI*sin(l*M_PI)/(2*dz)*fieldScale;
		Bfield[0]=Point[0]*Brtemp;
		Bfield[1]=Point[1]*Brtemp;
	} else {
		Bfield[0]=Bfield[1]=0.0;
	}
}
