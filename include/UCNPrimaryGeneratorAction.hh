#ifndef UCNPrimaryGeneratorAction_h
#define UCNPrimaryGeneratorAction_h 1

#include "UCNGlobals.hh"

#include "globals.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

#include "muParser.h"

class G4Event;
class G4ParticleGun;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class UCNPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:

  bool snapshotlog;
  bool tracklog;
  double trackloginterval;
  bool hitlog;
  bool flipspin;

  double xvar;

  int polarization;
  double t;
  double px, py, pz;
  double X ,Y, Z;

  double particleEnergy;
  double theta, phi;


  std::string sourcemode;
  double x_min, x_max, y_min, y_max, z_min, z_max;
  double r_min, r_max, phi_min, phi_max;//, z_min, z_max;
  double E_normal;
  double ActiveTime;
  bool PhaseSpaceWeighting;
  std::string sourcefile;

  std::string ParticleName;

  struct TParticleConfig{
    double tau; ///< lifetime
    double tmax; ///< max. simulation time
    double lmax; ///< max. trajectory length
    int polarization; ///< initial polarization
    double Emin; ///< min. initial energy
    double Emax; ///< max. initial energy
    mu::Parser spectrum; ///< Parsed energy spectrum given by user
    double phi_v_min; ///< Parsed minimum for initial azimuthal angle of velocity given by user
    double phi_v_max; ///< Parsed maximum for initial azimuthal angle of velocity given by user
    mu::Parser phi_v; ///< Parsed initial azimuthal angle distribution of velocity given by user
    double theta_v_min; ///< Parsed minimum for initial polarl angle of velocity given by user
    double theta_v_max; ///< Parsed maximum for initial polar angle of velocity given by user
    mu::Parser theta_v; ///< Parsed initial polar angle distribution of velocity given by user
  };

  std::map<std::string, TParticleConfig> pconfigs;

  TConfig geometryin;
  TConfig particlein;
  int ncount;
  //UCNPrimaryGeneratorAction(void);
  UCNPrimaryGeneratorAction(int simcount, int secondaries, TConfig GEOMIN, TConfig PARTIN);
  virtual ~UCNPrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event*);

private:

  G4ParticleGun* fParticleGun;
  //void ConfigInit(std::map<std::string, std::string> CONFIN);
  void   ConfigInit(TConfig CONFIN);
  void   TSource(TConfig CONFIN);

  double Spectrum(const std::string &particlename);
  void   AngularDist(const std::string &particlename, double &phi_v, double &theta_v);
  int    DicePolarisation(const std::string &particlename);
  void   RandomPointInSourceVolume();

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
