#ifndef UCNSteppingAction_H
#define UCNSteppingAction_H 1

#include "UCNTrackingAction.hh"
#include "UCNDetectorConstruction.hh"

#include "G4UserSteppingAction.hh"
#include "G4RunManager.hh"

#include <fstream>

class UCNSteppingAction : public G4UserSteppingAction
{
  public:

  UCNSteppingAction(int JOBNUM, std::string OUTPATH, int SECON, UCNTrackingAction* TAC, UCNDetectorConstruction* DTC, const bool *LOGINFO, double TRKLOGINT, const std::vector<double> &SNAPTIME);
    virtual ~UCNSteppingAction();

    virtual void UserSteppingAction(const G4Step*);

private:

  G4RunManager *runman;
  std::ofstream trackfile[3];
  std::ofstream hitfile[3];
  int pid;
  int secondaries;
  int jobnumber;
  std::string outpath;
  std::string name;
  int particle, polarisation;
  double t, x, y, z, vx, vy, vz, H, E;
  double B[4][4], Ei[3], V;
  double v1x, v1y, v1z, v2x, v2y, v2z, nx, ny, nz;
  int pol1, pol2, solid1, solid2;
  char phys_vol1[20], phys_vol2[20];
  std::string phys_name1, phys_name2;
  UCNDetectorConstruction* dtc;
  UCNTrackingAction* tac;
  bool Spinflip, hit;
  std::string pre_phys_name, post_phys_name;
  double pre_pol, post_pol;
  bool endlog, tracklog, hitlog, snapshotlog, spinlog;
  double trackloginterval;
  std::vector<double> snaptime;

  void OpenFile();
  void PrintData();
  void OpenHitFile();
  void PrintHitData();

};

#endif
