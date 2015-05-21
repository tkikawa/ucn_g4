#ifndef UCNActionInitialization_h
#define UCNActionInitialization_h 1

#include "UCNGlobals.hh"
#include "UCNDetectorConstruction.hh"

#include "G4VUserActionInitialization.hh"

/// Action initialization class.

class UCNActionInitialization : public G4VUserActionInitialization
{
public:
  int jobnumber;
  std::string outpath;

  TConfig particlein;
  TConfig geometryin;
  int secondaries;
  UCNDetectorConstruction* dtc;

  bool loginfo[5];
  double trackloginterval;
  std::vector<double> snaptime;


  UCNActionInitialization(int JOBNUM, std::string OUTPATH, int SECON, TConfig GEOMIN, TConfig PARIN, UCNDetectorConstruction* DTC);
  virtual ~UCNActionInitialization();
  
  virtual void BuildForMaster() const;
  virtual void Build() const;

  void ReadLogInfo(TConfig conf);
  
  //private:


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
