#include "ElectronBindingEnergy.hh"
#include "strutils.hh"
#include "SMExcept.hh"
#include <stdio.h>

const string BindingEnergyTable::shellnames = "KLMNOPQRST";

BindingEnergyTable::BindingEnergyTable(const Stringmap& m): Z(m.getDefault("Z",0)), nm(m.getDefault("name","")) {
	for(unsigned int n=0; n<shellnames.size(); n++) {
		vector<double> v;
		for(unsigned int s=1; s<=9; s++) {
			double b = m.getDefault(ctos(shellnames[n])+(n+s==1?"":itos(s)),0);
			if(b) v.push_back(b/1000.0);
			else break;
		}
		if(v.size()) eBinding.push_back(v);
		else break;
	}	
}

const vector<double>& BindingEnergyTable::getShellBinding(unsigned int n) const {
	if(n>=eBinding.size()) {
		SMExcept e("MissingShellInfo");
		e.insert("Z",Z);
		e.insert("name",nm);
		e.insert("shell",n);
		throw(e);
	}
	return eBinding[n];
}

double BindingEnergyTable::getSubshellBinding(unsigned int n, unsigned int m) const {
	const vector<double>& v = getShellBinding(n);
	if(m>=v.size()) {
		SMExcept e("MissingSubshellInfo");
		e.insert("Z",Z);
		e.insert("name",nm);
		e.insert("shell",n);
		e.insert("subshell",m);
		throw(e);
	}
	return v[m];
}

void BindingEnergyTable::display() const {
	printf("----- %i %s Electron Binding -----\n",Z,nm.c_str());
	for(unsigned int n=0; n<eBinding.size(); n++) {
		printf("\t%c:",shellnames[n]);
		for(unsigned int m=0; m<eBinding[n].size(); m++)
			printf("\t%.2f",eBinding[n][m]);
		printf("\n");
	}
}

//----------------------------------------------

BindingEnergyLibrary::BindingEnergyLibrary(const QFile& Q) {
	vector<Stringmap> v = Q.retrieve("binding");
	for(unsigned int i=0; i<v.size(); i++)
		tables.insert(std::pair<unsigned int,BindingEnergyTable*>(v[i].getDefault("Z",0),new BindingEnergyTable(v[i])));
}

BindingEnergyLibrary::~BindingEnergyLibrary() {
	for(map<unsigned int,BindingEnergyTable*>::const_iterator it = tables.begin(); it != tables.end(); it++)
		delete(it->second);
}

const BindingEnergyTable* BindingEnergyLibrary::getBindingTable(unsigned int Z) const {
	map<unsigned int,BindingEnergyTable*>::const_iterator it =  tables.find(Z);
	if(it==tables.end()) {
		SMExcept e("MissingElement");
		e.insert("Z",Z);
		throw(e);		
	}
	return it->second;
}

void BindingEnergyLibrary::display() const {
	for(map<unsigned int,BindingEnergyTable*>::const_iterator it = tables.begin(); it != tables.end(); it++)
		it->second->display();
}

