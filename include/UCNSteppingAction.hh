#ifndef UCNSteppingAction_H
#define UCNSteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

class UCNSteppingAction : public G4UserSteppingAction
{
  public:
    UCNSteppingAction();
    virtual ~UCNSteppingAction();

    virtual void UserSteppingAction(const G4Step*);
};

#endif

