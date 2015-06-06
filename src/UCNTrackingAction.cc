#include "G4SystemOfUnits.hh"

#include "UCNTrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4StepStatus.hh"
#include "G4VProcess.hh"
#include "G4ProcessType.hh"
#include "G4VPhysicalVolume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UCNTrackingAction::UCNTrackingAction(int JOBNUM, std::string OUTPATH, int SECON, UCNDetectorConstruction* DTC, bool ENDLOG)
{
  secondaries=SECON;
  jobnumber=JOBNUM;
  outpath = OUTPATH;
  dtc = DTC;
  runman =  G4RunManager::GetRunManager();
  issnapshot = false;
  endlog = ENDLOG;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UCNTrackingAction::~UCNTrackingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UCNTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  if(aTrack->GetParentID()!=0&&!secondaries){
    fpTrackingManager->SetStoreTrajectory(false);
    return;
  }
  else{
    fpTrackingManager->SetStoreTrajectory(true);
  }

  issnapshot = false;

  particle = runman->GetCurrentEvent()->GetEventID();
  tstart = aTrack->GetGlobalTime()/s;
  xstart = (aTrack->GetPosition()/m)[0];
  ystart = (aTrack->GetPosition()/m)[1];
  zstart = (aTrack->GetPosition()/m)[2];
  vxstart = (aTrack->GetVelocity()/(m/s))*(aTrack->GetMomentumDirection())[0];
  vystart = (aTrack->GetVelocity()/(m/s))*(aTrack->GetMomentumDirection())[1];
  vzstart = (aTrack->GetVelocity()/(m/s))*(aTrack->GetMomentumDirection())[2];
  if((aTrack->GetPolarization())[1]>0) polstart=1;
  else                                 polstart=-1;

  dtc->GetField()->GetCurrentFieldValue(tstart, xstart, ystart, zstart, B, Ei, V);
  Hstart = aTrack->GetKineticEnergy()/eV  + Epot(aTrack, V, polstart, B[3][0], zstart);
  Estart = aTrack->GetKineticEnergy()/eV;

  start_t = clock();
  NSpinflip = 0;
  Nhit = 0;
  Hmax = 0;
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

  dtc->GetField()->GetCurrentFieldValue(tend, xend, yend, zend, B, Ei, V);
  Hend = aTrack->GetKineticEnergy()/eV + Epot(aTrack, V, polend, B[3][0], zend);
  Eend = aTrack->GetKineticEnergy()/eV;

  stopID = 0;
  if(aTrack->GetStep()->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessType() == fDecay){
    stopID = -4;
  }
  else if(aTrack->GetStep()->GetPostStepPoint()->GetStepStatus() == fUserDefinedLimit){
    stopID = -1;
  }
  else if(aTrack->GetStep()->GetPostStepPoint()->GetStepStatus() == fWorldBoundary){
    stopID = -2;
  }
  else{
    //phys_name = aTrack->GetNextVolume()->GetName();
    phys_name = aTrack->GetStep()->GetPostStepPoint()->GetPhysicalVolume()->GetName();
    strcpy(phys_vol, phys_name.c_str());
    if(phys_vol != "World"){
      stopID = atoi(phys_vol+9);
    }
    else{
      stopID = 1;
    }
  }

  end_t = clock();
  ComputingTime = (double)(end_t - start_t)/CLOCKS_PER_SEC;

  Nstep      = aTrack->GetCurrentStepNumber();
  trajlength = aTrack->GetTrackLength()/m;

  spinflipprob=0;  

  pid = aTrack->GetDefinition()->GetPDGEncoding();
  if(pid==2112)      {pid=0;name="neutron";}
  else if(pid==2212) {pid=1;name="proton";}
  else if(pid==11)   {pid=2;name="electron";}
  else               return;

  if(endlog||issnapshot){
    if (!file[pid][issnapshot].is_open())OpenFile();
    PrintData();
  }

}

void UCNTrackingAction::SnapShotAction(const G4Track* aTrack)
{
  issnapshot=true;
  PostUserTrackingAction(aTrack);
  issnapshot=false;
}


void UCNTrackingAction::StepAction(double H, bool Spinflip, bool hit)
{
  if(Spinflip) NSpinflip++;
  if(hit) Nhit++;
  if(H>Hmax) Hmax = H;
}

void UCNTrackingAction::OpenFile(){
  std::string filesuffix;
  if(issnapshot) filesuffix = "snapshot.out";
  else           filesuffix = "end.out";
  std::ostringstream filename;
  filename << outpath << '/' << std::setw(12) << std::setfill('0') << jobnumber << name << filesuffix;
  std::cout << "Creating " << filename.str() << '\n';
  file[pid][issnapshot].open(filename.str().c_str());
  if (!file[pid][issnapshot].is_open()){
    std::cout << "Could not create" << filename.str() << '\n';
    exit(-1);
  }
  file[pid][issnapshot] << "jobnumber particle "
    "tstart xstart ystart zstart "
    "vxstart vystart vzstart "
    "polstart Hstart Estart "
    "tend xend yend zend "
    "vxend vyend vzend "
    "polend Hend Eend stopID Nspinflip spinflipprob "
    "ComputingTime Nhit Nstep trajlength Hmax\n";
  file[pid][issnapshot].precision(10);
}

void UCNTrackingAction::PrintData(){

  file[pid][issnapshot] << jobnumber << " " << particle << " "
			<< tstart << " " << xstart << " " << ystart << " " << zstart << " "
			<< vxstart << " " << vystart << " " << vzstart << " "
			<< polstart << " " << Hstart << " " << Estart << " "
			<< tend << " " << xend << " " << yend << " " << zend << " "
			<< vxend << " " << vyend << " " << vzend << " "
			<< polend << " " << Hend << " " << Eend << " " << stopID << " " << NSpinflip << " " << spinflipprob << " "
			<< ComputingTime << " " << Nhit << " " << Nstep << " " << trajlength << " " << Hmax << '\n';
  

}

double UCNTrackingAction::Epot(const G4Track* theTrack, double v, double pol, double b, double z){
  double epot = 
    (theTrack->GetDefinition()->GetPDGCharge()/coulomb)/ele_e*v
    - pol*(theTrack->GetDefinition()->GetPDGMagneticMoment()/(joule/tesla))/ele_e*b
    + (theTrack->GetDefinition()->GetPDGMass()/eV/c_0/c_0)*gravconst*z;
  return epot;
}
