#ifndef UCNTrackingAction_h
#define UCNTrackingAction_h 1

#include "G4UserTrackingAction.hh"


class UCNTrackingAction : public G4UserTrackingAction {

  public:
    UCNTrackingAction(){};
    virtual ~UCNTrackingAction(){};
   
    virtual void PreUserTrackingAction(const G4Track*);

};

#endif
