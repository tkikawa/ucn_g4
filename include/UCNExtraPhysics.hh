#ifndef UCNExtraPhysics_h
#define UCNExtraPhysics_h 1

#include "G4VPhysicsConstructor.hh"

class UCNExtraPhysics : public G4VPhysicsConstructor
{
public:

    UCNExtraPhysics();
    virtual ~UCNExtraPhysics();

    virtual void ConstructParticle();
    virtual void ConstructProcess();

private:

    void ConstructUCN();

    void AddBetaDecay();
};
#endif
