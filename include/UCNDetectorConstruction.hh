#ifndef UCNDetectorConstruction_h
#define UCNDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "G4UserLimits.hh"

#include "G4Field.hh"
#include "UCNGlobals.hh"
#include "UCN2DField.hh"
#include "UCN3DField.hh"
#include "UCNConductorField.hh"
#include "UCNField.hh"

// Open Asset Importer Library //
#include "assimp/Importer.hpp"
#include "assimp/scene.h"
#include "assimp/postprocess.h"

#include "G4TessellatedSolid.hh"

class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
#include "G4ThreeVector.hh"


#include "G4UCNMaterialPropertiesTable.hh"

class G4Material; 
class UCNField;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class UCNDetectorConstruction : public G4VUserDetectorConstruction
{
  public:

  G4UserLimits* g4limit;

  TConfig geometryin;

  std::string material_name;
  std::vector<std::string> materials;
  std::vector<TField*> fields;

  G4UCNMaterialPropertiesTable *mattbl[100];

  double a, cohcs, incohcs, scatcs, abscs, density, fermipot, loss;

  bool fieldIsInitialized;

  struct ignore{
    double ignorestart;
    double ignoreend;
    int stlid;
    int matid;
  };
  std::vector<ignore> ignores;

  UCNDetectorConstruction(double SIMTIME, TConfig GEOMIN);
  virtual ~UCNDetectorConstruction();

  public:
 
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

  private:

  G4VSolid *cad_solid[300], *cad_union[300];
  G4LogicalVolume *cad_logical[300];
  G4VPhysicalVolume *cad_physical[300];
 
  G4Material*     fVacuum;
  G4Material *ucn_material[100];

  Assimp::Importer importer;
  G4TessellatedSolid * volume_solid;
  const aiScene* scene;
  aiMesh* aim;
  G4VSolid* LoadCAD(char* filename);

public:
  
  UCNField* fField;
  double simtime;
  UCNField* GetField(){return fField;}
  double GetSimTime(){return simtime;};
  
private:
 
  void GetMaterial(int imat, std::string name);
  void GetElementValues(std::string name);
  void DefineMaterials();
  void ReadInField(TConfig conf);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
