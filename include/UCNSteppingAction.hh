#ifndef UCNSteppingAction_H
#define UCNSteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

#include <fstream>

class UCNSteppingAction : public G4UserSteppingAction
{
  public:

  std::ofstream trackfile;

  int jobnumber;

  int particle, polarisation;
  double t, x, y, z, vx, vy, vz, H, E;
  //double Bx, dBxdx, dBxdy, dBxdz, By, dBydx, dBydy, dBydz, Bz, dBzdx, dBzdy, dBzdz, Babs, dBdx, dBdy, dBdz, Ex, Ey, Ez, V;
  double B[4][4], Ei[3], V;


  UCNSteppingAction(int JOBNUM, std::string OUTPATH, std::string NAME);
    virtual ~UCNSteppingAction();

    virtual void UserSteppingAction(const G4Step*);
};

#endif

