#ifndef UCNPhysicsList_h
#define UCNPhysicsList_h 1

#include "G4VModularPhysicsList.hh"

class UCNPhysicsList: public G4VModularPhysicsList
{
public:

    UCNPhysicsList();
    virtual ~UCNPhysicsList();

    virtual void ConstructParticle();
    virtual void ConstructProcess();

    void SetCuts();

};

#endif
