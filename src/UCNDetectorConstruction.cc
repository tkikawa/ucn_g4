#include <iostream>
#include <sstream>
#include <string>

// USER //
#include "UCNDetectorConstruction.hh"

// CADMESH //
#include "CADMesh.hh"

// GEANT4 //
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4UCNMaterialPropertiesTable.hh"

// For field //
#include "G4Field.hh"
#include "UCNGlobals.hh"
#include "UCN2DField.hh"
#include "UCN3DField.hh"
#include "UCNConductorField.hh"
#include "UCNField.hh"

#include "G4Box.hh"
//#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UserLimits.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4UniformGravityField.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include "G4RepleteEofM.hh"
//#include "G4EqGravityField.hh"

#include "G4ClassicalRK4.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"
#include "G4PropagatorInField.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UCNDetectorConstruction::UCNDetectorConstruction(TConfig GEOMIN)
//UCNDetectorConstruction::UCNDetectorConstruction()
 : fVacuum(0)
{
  geometryin=GEOMIN;
  DefineMaterials();
  ReadInField(geometryin);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UCNDetectorConstruction::~UCNDetectorConstruction()
{
  if (fField) delete fField;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UCNDetectorConstruction::DefineMaterials()
{
  G4NistManager* nistMan = G4NistManager::Instance();

  fVacuum = nistMan->FindOrBuildMaterial("G4_Galactic");

  // --- Ni diffuse 10%
  /*
  G4UCNMaterialPropertiesTable* MPT = new G4UCNMaterialPropertiesTable();
  MPT->AddConstProperty("REFLECTIVITY",1.);
  MPT->AddConstProperty("DIFFUSION",0.1);
  MPT->AddConstProperty("FERMIPOT",252.0); // Gollub, Table 2.1 in neV
  MPT->AddConstProperty("SPINFLIP",0.);
  MPT->AddConstProperty("LOSS", 12.5e-5); //  Gollub, Table 2.1
  MPT->AddConstProperty("LOSSCS",0.);
  MPT->AddConstProperty("ABSCS",4.49); // 1/v loss cross-section  at room temp.
  MPT->AddConstProperty("SCATCS",18.5); // (incoherent) "elastic" scattering cs
  G4double neV = 1.e-9*eV;
  MPT->SetMicroRoughnessParameters(30*nm, 1*nm,
                                   180, 1000,
                                   0*degree, 90*degree,
                                   1*neV, 1000*neV,
                                   15, 15,
                                   0.01*degree);

  fGuideMaterial->SetMaterialPropertiesTable(MPT);
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  */

  int imat=0;
  double FermiReal; ///< Real part of Fermi potential
  double FermiImag; ///< Imaginary part of Fermi potential
  double DiffProb; ///< Diffuse reflection probability
  double SpinflipProb; ///< Probability for spin flip on reflection
  double RMSRoughness; ///< RMS roughness of surface, for MicroRoughness model reflections
  double CorrelLength; ///< Correlation length of surface roughness, for MicroRoughness model reflections
  bool UseMRModel; ///< Choose MicroRoughness model for reflections
  for (std::map<std::string, std::string>::iterator i = geometryin["MATERIALS"].begin(); i != geometryin["MATERIALS"].end(); i++){
    material_name = i->first;
    std::istringstream ss(i->second);
    ss >> FermiReal >> FermiImag >> DiffProb >> SpinflipProb >> RMSRoughness >> CorrelLength >> UseMRModel;

    if (ss){
      ucn_material[imat] = (G4Material*) nistMan->FindOrBuildMaterial(material_name);
      mattbl[imat] = new G4UCNMaterialPropertiesTable();
      mattbl[imat]->AddConstProperty("FERMIPOT",FermiReal);
      mattbl[imat]->AddConstProperty("LOSS",FermiImag/FermiReal);
      mattbl[imat]->AddConstProperty("DIFFUSION",DiffProb);
      mattbl[imat]->AddConstProperty("SPINFLIP",SpinflipProb);
      mattbl[imat]->AddConstProperty("MR_RRMS", RMSRoughness*nm);
      mattbl[imat]->AddConstProperty("MR_CORRLEN", CorrelLength*nm);
      //mattbl[imat]->AddConstProperty("",RMSRoughness);
      //mattbl[imat]->AddConstProperty("",CorrelLength);
      //mattbl[imat]->AddConstProperty("",UseMRModel);
      G4double neV = 1.e-9*eV;
      mattbl[imat]->SetMicroRoughnessParameters(CorrelLength*nm,
						RMSRoughness*nm,
						180, 1000,
						0*degree, 90*degree,
						1*neV, 1000*neV,
						15, 15,
						0.01*degree);
      ucn_material[imat]->SetMaterialPropertiesTable(mattbl[imat]);
      materials.push_back(material_name);
      imat++;

    }
    else
      std::cout << "Could not load material " << i->first << '\n';
  }





}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* UCNDetectorConstruction::Construct()
{
  //
  // World
  //

  G4double worldSizeX = 10.*m;
  G4double worldSizeY = 10.*m;
  G4double worldSizeZ = 10.*m;

  G4Box* solidWorld = new G4Box("World",
                                worldSizeX/2.,worldSizeY/2.,worldSizeZ/2.);

  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,
                                                    fVacuum,
                                                    "World");

  G4VPhysicalVolume* physiWorld = new G4PVPlacement(0,
                                                    G4ThreeVector(),
                                                    "World",
                                                    logicWorld,
                                                    0,
                                                    false,
                                                    0);

  logicWorld->SetVisAttributes(G4VisAttributes::Invisible);



  std::string line;
  std::string STLfile;
  std::string matname;
  int ID;
  int icad=0;
  char logical_name[100], physical_name[100];
  //char name[80];


  for (std::map<std::string, std::string>::iterator i = geometryin["GEOMETRY"].begin(); i != geometryin["GEOMETRY"].end(); i++){ // parse STLfile list
    //solid model;
    //std::istringstream(i->first) >> model.ID;
    std::istringstream(i->first) >> ID;
    //if (model.ID < 0){
    //std::cout << "You defined a solid with ID " << model.ID << " < 0! IDs have to be larger than zero!\n";
    if (ID < 0){
      std::cout << "You defined a solid with ID " << ID << " < 0! IDs have to be larger than zero!\n";
      exit(-1);
    }
    /*
    for (std::vector<solid>::iterator j = solids.begin(); j != solids.end(); j++){
      if (j->ID == model.ID){
	std::cout << "You defined solids with identical ID " << model.ID << "! IDs have to be unique!\n";
	exit(-1);
      }
    }
    */


    std::istringstream ss(i->second);
    ss >> STLfile >> matname;
    if (ss){
      for (unsigned j = 0; j < materials.size(); j++){
	if (matname == materials[j]){
	  //model.mat = materials[i];
	  //mat = materials[i];
	  //if (model.ID > 1){
	  if (ID > 1){
	    //mesh.ReadFile(STLfile.c_str(),solids.size(),name);
	    //model.name = name;
	    offset = G4ThreeVector(0, 0, 0);
	    std::cout<<"Reading: "<<STLfile<<std::endl;
	    CADMesh *mesh = new CADMesh((char*)STLfile.c_str(), "STL", mm, offset, false);
	    cad_solid[icad] = (G4VSolid*)mesh->TessellatedMesh();
	    sprintf(logical_name,"cad_logical_%d",icad);
	    cad_logical[icad] = new G4LogicalVolume(cad_solid[icad], ucn_material[j], logical_name, 0, 0, 0);
	    sprintf(physical_name,"cad_physical_%d",icad);
	    cad_physical[icad] = new G4PVPlacement(0, G4ThreeVector(), cad_logical[icad], physical_name, logicWorld, false, 0);
	    icad++;


	  }
	  //else
	  //model.name = "default solid";
	  break;
	}
	else if (j+1 == materials.size()){
	  printf("Material %s used for %s but not defined in geometry.in!",matname.c_str(), STLfile.c_str());
	  exit(-1);
	}
      }
    }
    else{
      // std::cout << "Could not load solid with ID " << model.ID << "! Did you define invalid parameters?\n";
      std::cout << "Could not load solid with ID " << ID << "! Did you define invalid parameters?\n";
      exit(-1);
    }


    double ignorestart, ignoreend;
    while (ss){
      ss >> ignorestart;
      if (!ss) // no more ignore times found                                                                                       
	break;
      char dash;
      ss >> dash;
      ss >> ignoreend;
      if (ss && dash == '-'){
	//model.ignoretimes.push_back(ignorestart);
	//model.ignoretimes.push_back(ignoreend);
      }
      else{
	std::cout << "Invalid ignoretimes in geometry.in" << '\n';
	exit(-1);
      }
    }

  }


  //
  //always return the physical World
  //
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal UCNField* UCNDetectorConstruction::fField = 0;


void UCNDetectorConstruction::ConstructSDandField()
{
  if (!fField) {

    fField = new UCNField(fields);
    fField->SetGravityActive(true);

     G4RepleteEofM* equation = new G4RepleteEofM(fField);
     // G4EqGravityField* equation = new G4EqGravityField(fField);

     G4FieldManager* fieldManager
      = G4TransportationManager::GetTransportationManager()->GetFieldManager();
     fieldManager->SetDetectorField(fField);


     // Set 12 to activate the spin tracking
     G4MagIntegratorStepper* stepper = new G4ClassicalRK4(equation,12);
     // G4MagIntegratorStepper* stepper = new G4ClassicalRK4(equation,8);

     G4double minStep           = 0.01*mm;

     G4ChordFinder* chordFinder = 
                   new G4ChordFinder((G4MagneticField*)fField,minStep,stepper);

     // Set accuracy parameters
     G4double deltaChord        = 3.0*mm;
     chordFinder->SetDeltaChord( deltaChord );

     G4double deltaOneStep      = 0.01*mm;
     fieldManager->SetAccuraciesWithDeltaOneStep(deltaOneStep);

     G4double deltaIntersection = 0.1*mm;
     fieldManager->SetDeltaIntersection(deltaIntersection);

     G4TransportationManager* transportManager =
                           G4TransportationManager::GetTransportationManager();

     G4PropagatorInField* fieldPropagator =
                                      transportManager->GetPropagatorInField();

     G4double epsMin            = 2.5e-7*mm;
     G4double epsMax            = 0.05*mm;

     fieldPropagator->SetMinimumEpsilonStep(epsMin);
     fieldPropagator->SetMaximumEpsilonStep(epsMax);

     fieldManager->SetChordFinder(chordFinder);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void UCNDetectorConstruction::ReadInField(TConfig conf){
  for (std::map<std::string, std::string>::iterator i = conf["FIELDS"].begin(); i != conf["FIELDS"].end(); i++){
    TField *f = NULL;
    std::istringstream ss(i->second);

    std::string ft;
    double Bscale, Escale, NullFieldTime, RampUpTime, FullFieldTime, RampDownTime;
    double Ibar, p1, p2, p3, p4, p5, p6;


    if (i->first == "2Dtable"){
      ss >> ft >> Bscale >> Escale >> NullFieldTime >> RampUpTime >> FullFieldTime >> RampDownTime;
      if (ss)f = new TabField(ft.c_str(), Bscale, Escale, NullFieldTime, RampUpTime, FullFieldTime, RampDownTime);
    }
    else if (i->first == "3Dtable"){
      double BoundaryWidth;
      ss >> ft >> Bscale >> Escale >> NullFieldTime >> RampUpTime >> FullFieldTime >> RampDownTime >> BoundaryWidth;
      if (ss)f = new TabField3(ft.c_str(), Bscale, Escale, NullFieldTime, RampUpTime, FullFieldTime, RampDownTime, BoundaryWidth);
    }
    else if (i->first == "InfiniteWireZ"){
      ss >> Ibar >> p1 >> p2;
      if (ss)f = new TInfiniteWireZ(p1, p2, Ibar);
    }
    else if (i->first == "InfiniteWireZCenter"){
      ss >> Ibar;
      if (ss)f = new TInfiniteWireZCenter(Ibar);
    }
    else if (i->first == "FiniteWire"){
      ss >> Ibar >> p1 >> p2 >> p3 >> p4 >> p5 >> p6;
      if (ss)f = new TFiniteWire(p1, p2, p3, p4, p5, p6, Ibar);
    }
    else if (i->first == "FiniteWireX"){
      ss >> Ibar >> p1 >> p2 >> p3;
      if (ss)f = new TFiniteWireX(p1, p2, p3, Ibar);
    }
    else if (i->first == "FiniteWireY"){
      ss >> Ibar >> p1 >> p2 >> p3;
      if (ss)f = new TFiniteWireY(p1, p2, p3, Ibar);
    }
    else if (i->first == "FiniteWireZ"){
      ss >> Ibar >> p1 >> p2 >> p3 >> p4;
      if (ss)f = new TFiniteWireZ(p1, p2, p3, p4, Ibar);
    }
    else if (i->first == "FiniteWireZCenter"){
      ss >> Ibar >> p1 >> p2;
      if (ss)f = new TFiniteWireZCenter(p1, p2, Ibar);
    }
    else if (i->first == "FullRacetrack"){
      ss >> Ibar >> p1 >> p2 >> p3;
      if (ss)f = new TFullRacetrack(p1, p2, p3, Ibar);
    }
    if (f)
      fields.push_back(f);
    else{
      std::cout << "Could not load field """ << i->first << """! Did you enter invalid parameters?\n";
      exit(-1);
    }
  }
}
