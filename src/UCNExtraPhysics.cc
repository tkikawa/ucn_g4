#include "UCNExtraPhysics.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4Transportation.hh"
#include "G4PhysicsListHelper.hh"

#include "G4UserSpecialCuts.hh"
#include "G4StepLimiter.hh"
#include "G4SystemOfUnits.hh"

#include "G4UCNLoss.hh"
#include "G4UCNAbsorption.hh"
#include "G4UCNMultiScattering.hh"
#include "G4UCNBoundaryProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UCNExtraPhysics::UCNExtraPhysics() 
    : G4VPhysicsConstructor("Extra") {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UCNExtraPhysics::~UCNExtraPhysics() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UCNExtraPhysics::ConstructParticle() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UCNExtraPhysics::ConstructProcess()
{
    aParticleIterator->reset();

    while ((*aParticleIterator)()) {
        G4ParticleDefinition* particle = aParticleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();

        if (!pmanager) {
            std::ostringstream o;
            o << "Particle " << particleName << "without a Process Manager";
            G4Exception("UCNExtraPhysics::ConstructProcess()","",
                         FatalException,o.str().c_str());
        }

        pmanager->AddDiscreteProcess(new G4StepLimiter());
        pmanager->AddDiscreteProcess(new G4UserSpecialCuts());
    }

    ConstructUCN();

    //activate spin tracking 
    G4Transportation::EnableUseMagneticMoment(true);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UCNExtraPhysics::ConstructUCN()
{
    aParticleIterator->reset();
    G4ProcessManager* pmanager = NULL;

    while ((*aParticleIterator)()) {
        G4ParticleDefinition* particle = aParticleIterator->value();
        pmanager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();

        if (!pmanager) {
	  std::ostringstream o;
	  o << "Particle " << particleName << "without a Process Manager";
	  G4Exception("UCNExtraPhysics::ConstructProcess()","",
		      FatalException,o.str().c_str());
        }

        if (particleName == "neutron") {
	  pmanager->AddDiscreteProcess(new G4UCNLoss());
	  pmanager->AddDiscreteProcess(new G4UCNAbsorption());
	  pmanager->AddDiscreteProcess(new G4UCNMultiScattering());
	  
	  G4UCNBoundaryProcess* ucnBoundaryProcess = 
	    new G4UCNBoundaryProcess();
	  ucnBoundaryProcess->SetMicroRoughness(true);
	  ucnBoundaryProcess->SetVerboseLevel(0);
	  
	  pmanager->AddDiscreteProcess(ucnBoundaryProcess);  
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
