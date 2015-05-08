#ifndef UCNRunAction_h
#define UCNRunAction_h 1

#include "G4Run.hh"

#include "G4UserRunAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class UCNRunAction : public G4UserRunAction
{
  public:
    UCNRunAction();
    virtual ~UCNRunAction();

  public:
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

  private:

};

#endif
