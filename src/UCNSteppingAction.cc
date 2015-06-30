#include "G4SystemOfUnits.hh"

#include "UCNSteppingAction.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

UCNSteppingAction::UCNSteppingAction(int JOBNUM, std::string OUTPATH, int SECON, UCNTrackingAction* TAC, UCNDetectorConstruction* DTC, const bool *LOGINFO, double TRKLOGINT, const std::vector<double> &SNAPTIME)
{
  secondaries=SECON;
  jobnumber=JOBNUM;
  outpath = OUTPATH;
  tac = TAC;
  dtc = DTC;
  endlog = LOGINFO[0];
  tracklog = LOGINFO[1];
  hitlog = LOGINFO[2];
  snapshotlog = LOGINFO[3];
  spinlog = LOGINFO[4];
  trackloginterval = TRKLOGINT;
  snaptime = SNAPTIME;

  runman =  G4RunManager::GetRunManager();
}

UCNSteppingAction::~UCNSteppingAction()
{;}

void UCNSteppingAction::UserSteppingAction(const G4Step * theStep)
{

  G4Track * theTrack = theStep->GetTrack();

  // check if it is alive
  if(theTrack->GetTrackStatus()!=fAlive) { return; }

  // check if it is primary
  if(theTrack->GetParentID()!=0&&!secondaries) { return; }

  particle = runman->GetCurrentEvent()->GetEventID();
  t = theTrack->GetGlobalTime()/s;
  x = (theTrack->GetPosition()/m)[0];
  y = (theTrack->GetPosition()/m)[1];
  z = (theTrack->GetPosition()/m)[2];
  vx = (theTrack->GetVelocity()/(m/s))*(theTrack->GetMomentumDirection())[0];
  vy = (theTrack->GetVelocity()/(m/s))*(theTrack->GetMomentumDirection())[1];
  vz = (theTrack->GetVelocity()/(m/s))*(theTrack->GetMomentumDirection())[2];
  if((theTrack->GetPolarization())[1]>0) polarisation=-1;
  else                                   polarisation=1;

  dtc->GetField()->GetCurrentFieldValue(t, x, y, z, B, Ei, V);
  
  H = theTrack->GetKineticEnergy()/eV + tac->Epot(theTrack, V, polarisation, B[3][0], x, y, z);
  E = theTrack->GetKineticEnergy()/eV;

  pre_pol = (theStep->GetPreStepPoint()->GetPolarization())[1];
  post_pol = (theStep->GetPostStepPoint()->GetPolarization())[1];
  pre_phys_name = theStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
  post_phys_name = theStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();
  if(pre_pol*post_pol<0) Spinflip = true;
  else Spinflip = false;
  if(pre_phys_name == "World" && post_phys_name != "World") hit = true;
  else hit = false;
  tac->StepAction(H, Spinflip, hit);

  pid = theTrack->GetDefinition()->GetPDGEncoding();
  if(pid==2112)      {pid=0;name="neutron";}
  else if(pid==2212) {pid=1;name="proton";}
  else if(pid==11)   {pid=2;name="electron";}
  else               return;

  if(hitlog){
    if(pre_phys_name != post_phys_name){
      
      v1x = (theStep->GetPreStepPoint()->GetVelocity()/(m/s))
	* (theStep->GetPreStepPoint()->GetMomentumDirection())[0];
      v1y = (theStep->GetPreStepPoint()->GetVelocity()/(m/s))
	* (theStep->GetPreStepPoint()->GetMomentumDirection())[1];
      v1z = (theStep->GetPreStepPoint()->GetVelocity()/(m/s))
	* (theStep->GetPreStepPoint()->GetMomentumDirection())[2];
      if((theStep->GetPreStepPoint()->GetPolarization())[1]>0) pol1 = -1;
      else pol1 = 1;

      v2x = (theStep->GetPostStepPoint()->GetVelocity()/(m/s))
	* (theStep->GetPostStepPoint()->GetMomentumDirection())[0];
      v2y = (theStep->GetPostStepPoint()->GetVelocity()/(m/s))
	* (theStep->GetPostStepPoint()->GetMomentumDirection())[1];
      v2z = (theStep->GetPostStepPoint()->GetVelocity()/(m/s))
	* (theStep->GetPostStepPoint()->GetMomentumDirection())[2];
      if((theStep->GetPostStepPoint()->GetPolarization())[1]>0) pol2 = -1;
      else pol2 = 1;

      phys_name1 = theStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
      strcpy(phys_vol1, phys_name1.c_str());
      solid1 = atoi(phys_vol1+9);

      phys_name2 = theStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();
      strcpy(phys_vol2, phys_name2.c_str());
      solid2 = atoi(phys_vol2+9);

      if (!hitfile[pid].is_open())OpenHitFile();
      PrintHitData();
    }
  }

  if(snapshotlog){
    double t0 = theStep->GetPreStepPoint()->GetGlobalTime()/s;
    double t1 = theStep->GetPostStepPoint()->GetGlobalTime()/s;
    for(unsigned int i=0;i<snaptime.size();i++){
      if(t0<snaptime[i]&&snaptime[i]<t1){
	tac->SnapShotAction(theTrack);
      }
    }
  }

  if(tracklog){
    if (!trackfile[pid].is_open())OpenFile();
    PrintData();
  }

}

void UCNSteppingAction::OpenFile(){

  std::string filesuffix = "track.out";
  std::ostringstream filename;
  filename << outpath << '/' << std::setw(12) << std::setfill('0') << jobnumber << name << filesuffix;
  std::cout << "Creating " << filename.str() << '\n';
  trackfile[pid].open(filename.str().c_str());
  if (!trackfile[pid].is_open()){
    std::cout << "Could not create" << filename.str() << '\n';
    exit(-1);
  }
  trackfile[pid] << "jobnumber particle polarisation "
    "t x y z vx vy vz "
    "H E Bx dBxdx dBxdy dBxdz By dBydx "
    "dBydy dBydz Bz dBzdx dBzdy dBzdz Babs dBdx dBdy dBdz Ex Ey Ez V\n";
  trackfile[pid].precision(10);
  
}

void UCNSteppingAction::PrintData(){

  trackfile[pid] << jobnumber << " " << particle << " " << polarisation << " "
	    << t << " " << x << " " << y << " " << z << " " << vx << " "
	    << vy << " " << vz << " " << H << " " << E << " ";
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      trackfile[pid] << B[i][j] << " ";
  trackfile[pid] << Ei[0] << " " << Ei[1] << " " << Ei[2] << " " << V << '\n';

}


void UCNSteppingAction::OpenHitFile(){

  std::string filesuffix = "hit.out";
  std::ostringstream filename;
  filename << outpath << '/' << std::setw(12) << std::setfill('0') << jobnumber << name << filesuffix;
  std::cout << "Creating " << filename.str() << '\n';
  hitfile[pid].open(filename.str().c_str());
  if (!hitfile[pid].is_open()){
    std::cout << "Could not create" << filename.str() << '\n';
    exit(-1);
  }
  hitfile[pid] << "jobnumber particle "
    "t x y z v1x v1y v1z pol1 "
    "v2x v2y v2z pol2 "
    "nx ny nz solid1 solid2\n";
  hitfile[pid].precision(10);
  
}

void UCNSteppingAction::PrintHitData(){

  hitfile[pid] << jobnumber << " " << particle << " " << polarisation << " "
	       << t << " " << x << " " << y << " " << z << " " << v1x << " "
	       << v1y << " " << v1z << " " << pol1 << " "<< v1x << " "
	       << v2y << " " << v2z << " " << pol2 << " "<<nx<< " "<<ny<< " "
	       <<nz<< " " <<solid1<< " " <<solid2<< '\n';

}
