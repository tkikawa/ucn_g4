#ifndef UCNTrackingAction_h
#define UCNTrackingAction_h 1

#include "G4UserTrackingAction.hh"

#include <fstream>

class UCNTrackingAction : public G4UserTrackingAction {

  public:
  std::ofstream file;

  int jobnumber;
  double tstart, xstart, ystart, zstart, vxstart, vystart, vzstart, Hstart, Estart, tend, xend, yend, zend, vxend, vyend, vzend, Hend, Eend, spinflipprob, ComputingTime, trajlength, Hmax;
  int particle, polstart, polend, stopID, NSpinflip, Nhit, Nstep;


  UCNTrackingAction(int JOBNUM, std::string OUTPATH, std::string NAME);
  virtual ~UCNTrackingAction();
   
  virtual void PreUserTrackingAction(const G4Track*);
  virtual void PostUserTrackingAction(const G4Track*);

};

#endif
