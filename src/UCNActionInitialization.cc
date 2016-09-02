#include "UCNActionInitialization.hh"
#include "UCNPrimaryGeneratorAction.hh"
#include "UCNRunAction.hh"
#include "UCNSteppingAction.hh"
#include "UCNTrackingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UCNActionInitialization::UCNActionInitialization(unsigned long long int JOBNUM, std::string OUTPATH, int SECON, TConfig GEOMIN, TConfig PARTIN, UCNDetectorConstruction* DTC)
 : G4VUserActionInitialization()
{
  jobnumber=JOBNUM;
  outpath=OUTPATH;
  secondaries=SECON;
  particlein=PARTIN;
  geometryin=GEOMIN;
  dtc=DTC;
  ReadLogInfo(particlein);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UCNActionInitialization::~UCNActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UCNActionInitialization::BuildForMaster() const
{
  UCNRunAction* rac = new UCNRunAction();
  SetUserAction(rac);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UCNActionInitialization::Build() const
{
  UCNPrimaryGeneratorAction* gac = new UCNPrimaryGeneratorAction(geometryin, particlein);
  SetUserAction(gac);
  UCNRunAction* rac = new UCNRunAction();
  SetUserAction(rac);
  UCNTrackingAction* tac = new UCNTrackingAction(jobnumber, outpath, secondaries, dtc, loginfo[0]);
  SetUserAction(tac);
  UCNSteppingAction* sac = new UCNSteppingAction(jobnumber, outpath, secondaries, tac, dtc, loginfo, trackloginterval, snaptime);
  SetUserAction(sac);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UCNActionInitialization::ReadLogInfo(TConfig conf)
{
  bool endlog = false;
  bool tracklog = false;
  bool hitlog = false;
  bool snapshotlog = false;
  bool spinlog = false;
  trackloginterval = 1e-3;

  double sstime = -1;
  std::istringstream(conf["all"]["endlog"]) >> endlog;
  std::istringstream(conf["all"]["tracklog"]) >> tracklog;
  std::istringstream(conf["all"]["trackloginterval"]) >> trackloginterval;
  std::istringstream(conf["all"]["hitlog"]) >> hitlog;
  std::istringstream(conf["all"]["snapshotlog"]) >> snapshotlog;
  std::istringstream snapshots(conf["all"]["snapshots"]);
  if (snapshotlog){
    while(snapshots >> sstime){snaptime.push_back(sstime);}
  }
  std::istringstream(conf["all"]["spintlog"]) >> spinlog;
  
  loginfo[0] = endlog;
  loginfo[1] = tracklog;
  loginfo[2] = hitlog;
  loginfo[3] = snapshotlog;
  loginfo[4] = spinlog;
}
