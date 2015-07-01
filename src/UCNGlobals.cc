#include <iostream>
#include <fstream>
#include <cmath>

#include "UCNGlobals.hh"

// print progress in percent
void PrintPercent(double percentage, int &lastprint){
	// write status to console
	// one point per 2 percent of endtime
	while (lastprint < percentage*100){
		lastprint += 2;
		if (lastprint % 10 == 0)
			std::cout << lastprint << "%";
		else
			std::cout << ".";
//		std::cout.flush();
	}
}

// energy distribution of protons (0 < 'Energie' < 750 eV)
// proton recoil spectrum from "Diplomarbeit M. Simson"
// result always < 1!
double ProtonBetaSpectrum(double E){
	double DeltaM = m_n - m_p;
	double Xi = m_e / DeltaM;
	double Sigma = 1 - 2*E*m_n / pow(DeltaM, 2) / pow(c_0, 2);
	double g1 = pow(((Sigma - pow(Xi, 2)) / Sigma), 2) * sqrt(1 - Sigma) * (4*(1 + pow(Xi, 2) / Sigma) - 				4/3*((Sigma - pow(Xi, 2)) / Sigma * (1 - Sigma)));
	double g2 = pow(((Sigma - pow(Xi, 2)) / Sigma), 2) * sqrt(1 - Sigma) * (4*(1 + pow(Xi, 2) / Sigma - 2 * Sigma) - 	4/3*((Sigma - pow(Xi, 2)) / Sigma * (1 - Sigma)));
	double littlea = -0.1017;
	return 0.5*(g1 + littlea * g2);
}


// energy distribution of electrons (0 < 'Energie' < 782 keV)
// from "http://hyperphysics.phy-astr.gsu.edu/Hbase/nuclear/beta2.html"
// result always < 1!
double ElectronBetaSpectrum(double E){
	double Qvalue = 0.782; //[MeV]
	return 8.2*sqrt(E*1e-6*E*1e-6 + 2*E*1e-6*m_e*c_0*c_0*1e-6) * pow(Qvalue - E*1e-6, 2) * (E*1e-6 + m_e*c_0*c_0*1e-6);
}


typedef std::map<std::string, std::map<std::string, std::string> > TConfig;


//read variables from *.in file into map
void ReadInFile(const char *inpath, TConfig &vars){
	std::ifstream infile(inpath);
	char c;
	std::string rest,section,key;
	while (infile && (infile >> std::ws) && (c = infile.peek())){
		if (c == '[' && infile.ignore()){
			if (infile.peek() == '/'){
				section = "";
			}
			else{
				std::getline(infile, section, ']');
//				std::cout << "\nsection: " << section.c_str() << '\n';
			}
			std::getline(infile,rest);
		}
		else if (c == '#')
			std::getline(infile,rest);
		else if (section != ""){
			infile >> key;
			std::getline(infile,rest);
			if (infile){
				std::string::size_type l = rest.find('#');
				if (l == std::string::npos)
					vars[section][key] = rest;
				else
					vars[section][key] = rest.substr(0,l);
//				std::cout << key << " " << vars[section][key] << '\n';
			}
		}
		else
			std::getline(infile,rest);
	}
}
