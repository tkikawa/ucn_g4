#include "UCNRunAction.hh"

#include "G4Run.hh"

#include "G4ProcessTable.hh"
#include "G4UCNBoundaryProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UCNRunAction::UCNRunAction() : G4UserRunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UCNRunAction::~UCNRunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UCNRunAction::BeginOfRunAction(const G4Run*)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UCNRunAction::EndOfRunAction(const G4Run*)
{
  G4ProcessTable* processTable = G4ProcessTable::GetProcessTable();

  G4UCNBoundaryProcess* ucnBoundaryProcess = 
             (G4UCNBoundaryProcess*) processTable->
             FindProcess("UCNBoundaryProcess",G4Neutron::NeutronDefinition());

  ucnBoundaryProcess->BoundaryProcessSummary();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
