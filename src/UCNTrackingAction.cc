#include "G4SystemOfUnits.hh"

#include "UCNTrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UCNTrackingAction::UCNTrackingAction(int JOBNUM, std::string OUTPATH, int SECON)
{
  secondaries=SECON;
  jobnumber=JOBNUM;
  outpath = OUTPATH;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UCNTrackingAction::~UCNTrackingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UCNTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  // Create trajectory only for primaries
  //if(aTrack->GetParentID()==0)
  //{ fpTrackingManager->SetStoreTrajectory(true); }
  //else
  //{ fpTrackingManager->SetStoreTrajectory(false); }

  if(aTrack->GetParentID()!=0&&!secondaries){
    fpTrackingManager->SetStoreTrajectory(false);
    return;
  }
  else{
    fpTrackingManager->SetStoreTrajectory(true);
  }

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

  if(aTrack->GetParentID()!=0&&!secondaries) { return; }

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

  pid = aTrack->GetDefinition()->GetPDGEncoding();
  if(pid==2112)      {pid=0;name="neutron";}
  else if(pid==2212) {pid=1;name="proton";}
  else if(pid==11)   {pid=2;name="electron";}
  else               return;

  if (!file[pid].is_open())OpenFile();
  PrintData();

}

void UCNTrackingAction::OpenFile(){

  std::string filesuffix = "end.out";
  std::ostringstream filename;
  filename << outpath << '/' << std::setw(12) << std::setfill('0') << jobnumber << name << filesuffix;
  std::cout << "Creating " << filename.str() << '\n';
  file[pid].open(filename.str().c_str());
  if (!file[pid].is_open()){
    std::cout << "Could not create" << filename.str() << '\n';
    exit(-1);
  }
  file[pid] << "jobnumber particle "
    "tstart xstart ystart zstart "
    "vxstart vystart vzstart "
    "polstart Hstart Estart "
    "tend xend yend zend "
    "vxend vyend vzend "
    "polend Hend Eend stopID Nspinflip spinflipprob "
    "ComputingTime Nhit Nstep trajlength Hmax\n";
  file[pid].precision(10);
}

void UCNTrackingAction::PrintData(){

  file[pid] << jobnumber << " " << particle << " "
	    << tstart << " " << xstart << " " << ystart << " " << zstart << " "
	    << vxstart << " " << vystart << " " << vzstart << " "
	    << polstart << " " << Hstart << " " << Estart << " "
	    << tend << " " << xend << " " << yend << " " << zend << " "
	    << vxend << " " << vyend << " " << vzend << " "
	    << polend << " " << Hend << " " << Eend << " " << stopID << " " << NSpinflip << " " << spinflipprob << " "
	    << ComputingTime << " " << Nhit << " " << Nstep << " " << trajlength << " " << Hmax << '\n';
  

}
