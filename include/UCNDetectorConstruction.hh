#ifndef UCNDetectorConstruction_h
#define UCNDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "G4Field.hh"
#include "UCNGlobals.hh"
#include "UCN2DField.hh"
#include "UCN3DField.hh"
#include "UCNConductorField.hh"
#include "UCNField.hh"

class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
#include "G4ThreeVector.hh"


#include "G4UCNMaterialPropertiesTable.hh"

class G4Material; 
//class G4UniformGravityField;
class UCNField;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class UCNDetectorConstruction : public G4VUserDetectorConstruction
{
  public:

  TConfig geometryin;

  std::string material_name;
  std::vector<std::string> materials;

  std::vector<TField*> fields;

  G4UCNMaterialPropertiesTable *mattbl[100];


  UCNDetectorConstruction(TConfig GEOMIN);
  //UCNDetectorConstruction();
  virtual ~UCNDetectorConstruction();

  public:
 
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

  private:


  G4ThreeVector offset;
  G4VSolid *cad_solid[300];
  G4LogicalVolume *cad_logical[300];
  G4VPhysicalVolume *cad_physical[300];
 
  G4Material*     fVacuum;
  G4Material *ucn_material[100];
  
  static G4ThreadLocal UCNField* fField;

  
private:
 
  void GetMaterial(int imat, std::string name);
  void DefineMaterials();
  void ReadInField(TConfig conf);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
