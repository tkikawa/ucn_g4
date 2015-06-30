#ifndef UCNTrackingAction_h
#define UCNTrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "G4RunManager.hh"

#include "UCNDetectorConstruction.hh"

#include <fstream>
#include <sys/time.h>

class UCNTrackingAction : public G4UserTrackingAction {

  public:
  
  UCNTrackingAction(int JOBNUM, std::string OUTPATH, int SECON, UCNDetectorConstruction* DTC, bool ENDLOG);
  virtual ~UCNTrackingAction();
   
  virtual void PreUserTrackingAction(const G4Track*);
  virtual void PostUserTrackingAction(const G4Track*);

  void SnapShotAction(const G4Track*);
  void StepAction(double H, bool Spinflip, bool hit);
  double Epot(const G4Track* theTrack, double v, double pol, double b, double x, double y, double z);
  int SolidID(G4String phys_name);

private:

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
  double tstart, xstart, ystart, zstart, vxstart, vystart, vzstart, Hstart, Estart, Bstart, Ustart, tend, xend, yend, zend, vxend, vyend, vzend, Hend, Eend, Bend, Uend, spinflipprob, ComputingTime, trajlength, Hmax;
  int particle, polstart, polend, stopID, NSpinflip, Nhit, Nstep, solidstart, solidend;

  double B[4][4], Ei[3], V;
  clock_t start_t, end_t;
  char phys_vol[20];
  UCNDetectorConstruction* dtc;

  void OpenFile();
  void PrintData();
};

#endif
