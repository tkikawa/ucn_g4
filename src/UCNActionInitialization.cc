#include "UCNActionInitialization.hh"
#include "UCNPrimaryGeneratorAction.hh"
#include "UCNRunAction.hh"
#include "UCNSteppingAction.hh"
#include "UCNTrackingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UCNActionInitialization::UCNActionInitialization(int JOBNUM, std::string OUTPATH, int SECON, TConfig GEOMIN, TConfig PARTIN)
 : G4VUserActionInitialization()
{
  jobnumber=JOBNUM;
  outpath=OUTPATH;
  secondaries=SECON;
  particlein=PARTIN;
  geometryin=GEOMIN;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UCNActionInitialization::~UCNActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UCNActionInitialization::BuildForMaster() const
{
  SetUserAction(new UCNRunAction());
  SetUserAction(new UCNSteppingAction(jobnumber, outpath, secondaries));
  SetUserAction(new UCNTrackingAction(jobnumber, outpath, secondaries));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UCNActionInitialization::Build() const
{
  SetUserAction(new UCNPrimaryGeneratorAction(secondaries, geometryin, particlein));
  //SetUserAction(new UCNPrimaryGeneratorAction());

  SetUserAction(new UCNRunAction());
  SetUserAction(new UCNSteppingAction(jobnumber, outpath, secondaries));
  SetUserAction(new UCNTrackingAction(jobnumber, outpath, secondaries));
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
