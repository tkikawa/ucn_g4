#ifndef UCNTrackingAction_h
#define UCNTrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "G4RunManager.hh"

#include "UCNDetectorConstruction.hh"

#include <fstream>
#include <sys/time.h>

class UCNTrackingAction : public G4UserTrackingAction {

  public:
  
  G4RunManager *runman;

  std::ofstream file[3][2];
  int pid;
  bool issnapshot;

  bool endlog;
  
  int secondaries;
  int jobnumber;
  std::string outpath;
  std::string name;

  double simtime;
  
  double tstart, xstart, ystart, zstart, vxstart, vystart, vzstart, Hstart, Estart, tend, xend, yend, zend, vxend, vyend, vzend, Hend, Eend, spinflipprob, ComputingTime, trajlength, Hmax;
  int particle, polstart, polend, stopID, NSpinflip, Nhit, Nstep;

  double B[4][4], Ei[3], V;

  clock_t start_t, end_t;

  char phys_vol[20];
  std::string phys_name;

  UCNDetectorConstruction* dtc;

  UCNTrackingAction(int JOBNUM, std::string OUTPATH, int SECON, UCNDetectorConstruction* DTC, bool ENDLOG);
  virtual ~UCNTrackingAction();
   
  virtual void PreUserTrackingAction(const G4Track*);
  virtual void PostUserTrackingAction(const G4Track*);

  void SnapShotAction(const G4Track*);

  void StepAction(double H, bool Spinflip, bool hit);
  void OpenFile();
  void PrintData();
  double Epot(const G4Track* theTrack, double v, double pol, double b, double z);
};

#endif
