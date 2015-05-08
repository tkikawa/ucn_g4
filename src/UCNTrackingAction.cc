#include "G4SystemOfUnits.hh"

#include "UCNTrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UCNTrackingAction::UCNTrackingAction(int JOBNUM, std::string OUTPATH, std::string NAME)
{
  jobnumber=JOBNUM;
  std::string filesuffix = "end.out";
  std::string outpath = OUTPATH;
  std::string name = NAME;

  if (!file.is_open()){
    std::ostringstream filename;
    filename << outpath << '/' << std::setw(12) << std::setfill('0') << jobnumber << name << filesuffix;
    std::cout << "Creating " << filename.str() << '\n';
    file.open(filename.str().c_str());
    if (!file.is_open()){
      std::cout << "Could not create" << filename.str() << '\n';
      exit(-1);
    }
    file << "jobnumber particle "
      "tstart xstart ystart zstart "
      "vxstart vystart vzstart "
      "polstart Hstart Estart "
      "tend xend yend zend "
      "vxend vyend vzend "
      "polend Hend Eend stopID Nspinflip spinflipprob "
      "ComputingTime Nhit Nstep trajlength Hmax\n";
    file.precision(10);
  }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UCNTrackingAction::~UCNTrackingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UCNTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  // Create trajectory only for primaries
  if(aTrack->GetParentID()==0)
  { fpTrackingManager->SetStoreTrajectory(true); }
  else
  { fpTrackingManager->SetStoreTrajectory(false); }

  particle=0;
  tstart = aTrack->GetGlobalTime()/s;
  xstart = (aTrack->GetPosition()/m)[0];
  ystart = (aTrack->GetPosition()/m)[1];
  zstart = (aTrack->GetPosition()/m)[2];
  vxstart = (aTrack->GetVelocity()/(m/s))*(aTrack->GetMomentumDirection())[0];
  vystart = (aTrack->GetVelocity()/(m/s))*(aTrack->GetMomentumDirection())[1];
  vzstart = (aTrack->GetVelocity()/(m/s))*(aTrack->GetMomentumDirection())[2];
  if((aTrack->GetPolarization())[1]>0) polstart=1;
  else                                 polstart=-1;
  Hstart = aTrack->GetTotalEnergy()/eV;
  Estart = aTrack->GetKineticEnergy()/eV;



}


void UCNTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  tend = aTrack->GetGlobalTime()/s;
  xend = (aTrack->GetPosition()/m)[0];
  yend = (aTrack->GetPosition()/m)[1];
  zend = (aTrack->GetPosition()/m)[2];
  vxend = (aTrack->GetVelocity()/(m/s))*(aTrack->GetMomentumDirection())[0];
  vyend = (aTrack->GetVelocity()/(m/s))*(aTrack->GetMomentumDirection())[1];
  vzend = (aTrack->GetVelocity()/(m/s))*(aTrack->GetMomentumDirection())[2];
  if((aTrack->GetPolarization())[1]>0) polend=1;
  else                                 polend=-1;
  Hend = aTrack->GetTotalEnergy()/eV;
  Eend = aTrack->GetKineticEnergy()/eV;

  stopID=0;
  NSpinflip=0;
  spinflipprob=0;
  ComputingTime=0;
  Nhit=0;

  Nstep      = aTrack->GetCurrentStepNumber();
  trajlength = aTrack->GetTrackLength()/m;

  Hmax=0;


  file    << jobnumber << " " << particle << " "
	  << tstart << " " << xstart << " " << ystart << " " << zstart << " "
	  << vxstart << " " << vystart << " " << vzstart << " "
	  << polstart << " " << Hstart << " " << Estart << " "
	  << tend << " " << xend << " " << yend << " " << zend << " "
	  << vxend << " " << vyend << " " << vzend << " "
	  << polend << " " << Hend << " " << Eend << " " << stopID << " " << NSpinflip << " " << spinflipprob << " "
	  << ComputingTime << " " << Nhit << " " << Nstep << " " << trajlength << " " << Hmax << '\n';


}
