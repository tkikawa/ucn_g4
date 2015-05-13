#include "UCNActionInitialization.hh"
#include "UCNPrimaryGeneratorAction.hh"
#include "UCNRunAction.hh"
#include "UCNSteppingAction.hh"
#include "UCNTrackingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UCNActionInitialization::UCNActionInitialization(int JOBNUM, std::string OUTPATH, TConfig CONFIN, TConfig GEOMIN, TConfig PARTIN)
 : G4VUserActionInitialization()
{
  jobnumber=JOBNUM;
  outpath=OUTPATH;

  simtype=1;
  ConfigInit(CONFIN);
  if(simtype!=1){
    std::cout<<"Only the particle simualtion is available.\n";
    exit(-1);
  }
  particlein=PARTIN;
  geometryin=GEOMIN;

  std::istringstream sourceconf(geometryin["SOURCE"].begin()->second);
  sourceconf >> ParticleName;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UCNActionInitialization::~UCNActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UCNActionInitialization::BuildForMaster() const
{
  SetUserAction(new UCNRunAction());
  SetUserAction(new UCNSteppingAction(jobnumber, outpath, secondaries, ParticleName));
  SetUserAction(new UCNTrackingAction(jobnumber, outpath, secondaries, ParticleName));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UCNActionInitialization::Build() const
{
  SetUserAction(new UCNPrimaryGeneratorAction(simcount, secondaries, geometryin, particlein));
  //SetUserAction(new UCNPrimaryGeneratorAction());

  SetUserAction(new UCNRunAction());
  SetUserAction(new UCNSteppingAction(jobnumber, outpath, secondaries, ParticleName));
  SetUserAction(new UCNTrackingAction(jobnumber, outpath, secondaries, ParticleName));
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void UCNActionInitialization::ConfigInit(TConfig config){
  /* setting default values */
  simtype = PARTICLE;
  neutdist = 0;
  simcount = 1;
  /*end default values*/

  /* read variables from map by casting strings in map into istringstreams and extracting value with ">>"-operator */
  std::istringstream(config["global"]["simtype"])>> simtype;
  std::istringstream(config["global"]["neutdist"])>> neutdist;
  

  std::istringstream(config["global"]["simcount"])>> simcount;
  std::istringstream(config["global"]["simtime"])>> SimTime;
  std::istringstream(config["global"]["secondaries"])>> secondaries;
  std::istringstream(config["global"]["BCutPlane"])>> BCutPlanePoint[0] >> BCutPlanePoint[1] >> BCutPlanePoint[2]
					      >> BCutPlanePoint[3] >> BCutPlanePoint[4] >> BCutPlanePoint[5]
					      >> BCutPlanePoint[6] >> BCutPlanePoint[7] >> BCutPlanePoint[8]
					      >> BCutPlaneSampleCount1 >> BCutPlaneSampleCount2;
}
