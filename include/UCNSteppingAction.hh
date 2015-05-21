#ifndef UCNSteppingAction_H
#define UCNSteppingAction_H 1

#include "G4UserSteppingAction.hh"
#include "G4RunManager.hh"

#include "UCNTrackingAction.hh"
#include "UCNDetectorConstruction.hh"

#include <fstream>

class UCNSteppingAction : public G4UserSteppingAction
{
  public:

  G4RunManager *runman;

  std::ofstream trackfile[3];
  int pid;

  int secondaries;
  int jobnumber;
  std::string outpath;
  std::string name;

  int particle, polarisation;
  double t, x, y, z, vx, vy, vz, H, E;
  double B[4][4], Ei[3], V;

  UCNDetectorConstruction* dtc;

  UCNTrackingAction* tac;
  bool Spinflip, hit;
  std::string pre_phys_name, post_phys_name;
  double pre_pol, post_pol;

  bool endlog, tracklog, hitlog, snapshotlog, spinlog;
  double trackloginterval;
  std::vector<double> snaptime;

  UCNSteppingAction(int JOBNUM, std::string OUTPATH, int SECON, UCNTrackingAction* TAC, UCNDetectorConstruction* DTC, const bool *LOGINFO, double TRKLOGINT, const std::vector<double> &SNAPTIME);
    virtual ~UCNSteppingAction();

    virtual void UserSteppingAction(const G4Step*);

  void OpenFile();
  void PrintData();

};

#endif
