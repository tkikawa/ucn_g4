#ifndef UCNActionInitialization_h
#define UCNActionInitialization_h 1

#include "UCNGlobals.hh"

#include "G4VUserActionInitialization.hh"

/// Action initialization class.

class UCNActionInitialization : public G4VUserActionInitialization
{
public:
  int jobnumber;
  std::string outpath;

  TConfig particlein;
  TConfig geometryin;

  std::string ParticleName;

  double SimTime; ///< max. simulation time
  int neutdist;
  int simcount; ///< number of particles for MC simulation (read from config)
  int simtype ; ///< type of particle which shall be simulated (read from config)
  int secondaries ; ///< should secondary particles be simulated? (read from config)
  double BCutPlanePoint[9]; ///< 3 points on plane for field slice (read from config)
  int BCutPlaneSampleCount1; ///< number of field samples in BCutPlanePoint[3..5]-BCutPlanePoint[0..2] direction (read from config)
  int BCutPlaneSampleCount2; ///< number of field samples in BCutPlanePoint[6..8]-BCutPlanePoint[0..2] direction (read from config)



  UCNActionInitialization(int JOBNUM, std::string OUTPATH, TConfig CONFMIN, TConfig GEOMIN, TConfig PARIN);
  virtual ~UCNActionInitialization();
  
  virtual void BuildForMaster() const;
  virtual void Build() const;
  
private:
  void ConfigInit(TConfig CONFIN);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
