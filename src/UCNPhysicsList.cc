#include "UCNPhysicsList.hh"

#include "G4ParticleTypes.hh"

#include "G4DecayPhysics.hh"
#include "UCNExtraPhysics.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UCNPhysicsList::UCNPhysicsList() : G4VModularPhysicsList() 
{
    RegisterPhysics(new G4DecayPhysics());
    RegisterPhysics(new UCNExtraPhysics());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UCNPhysicsList::~UCNPhysicsList() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UCNPhysicsList::ConstructParticle()
{
    G4Neutron::NeutronDefinition();
    G4Proton::ProtonDefinition();
    G4Electron::ElectronDefinition();
    G4AntiNeutrinoE::AntiNeutrinoEDefinition();
    G4MuonPlus::MuonPlusDefinition();
    G4MuonMinus::MuonMinusDefinition();

    G4GenericIon::GenericIonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UCNPhysicsList::ConstructProcess()
{
    G4VModularPhysicsList::ConstructProcess();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UCNPhysicsList::SetCuts()
{
    SetCutsWithDefault();
}
