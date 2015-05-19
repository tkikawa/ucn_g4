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
#include "G4UnionSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

//#include "G4GeometryManager.hh"
//#include "G4PhysicalVolumeStore.hh"
//#include "G4LogicalVolumeStore.hh"
//#include "G4SolidStore.hh"

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

UCNDetectorConstruction::UCNDetectorConstruction(int SIMTIME, TConfig GEOMIN)
//UCNDetectorConstruction::UCNDetectorConstruction()
 : fVacuum(0)
{
  geometryin=GEOMIN;
  DefineMaterials();
  ReadInField(geometryin);
  g4limit = new G4UserLimits(DBL_MAX,DBL_MAX,SIMTIME*s);
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

  G4double neV = 1.e-9*eV;

  int imat = 0;
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

    if (material_name=="default"||material_name=="Default")continue;
    if (ss){
      std::cout<<"Material setting for "<<material_name.c_str()<<"..."<<std::flush;
      GetMaterial(imat, material_name);

      mattbl[imat] = new G4UCNMaterialPropertiesTable();
      mattbl[imat]->AddConstProperty("FERMIPOT",FermiReal);
      mattbl[imat]->AddConstProperty("LOSS",FermiImag/FermiReal);
      mattbl[imat]->AddConstProperty("DIFFUSION",DiffProb);
      mattbl[imat]->AddConstProperty("SPINFLIP",SpinflipProb);
      mattbl[imat]->AddConstProperty("MR_RRMS", RMSRoughness*m);
      mattbl[imat]->AddConstProperty("MR_CORRLEN", CorrelLength*m);

      mattbl[imat]->AddConstProperty("REFLECTIVITY",1.);
      mattbl[imat]->AddConstProperty("LOSSCS",0.);
      mattbl[imat]->AddConstProperty("ABSCS",abscs);
      mattbl[imat]->AddConstProperty("SCATCS",scatcs);

      mattbl[imat]->SetMicroRoughnessParameters(CorrelLength*m,
						RMSRoughness*m,
						180, 1000,
						0*degree, 90*degree,
						1*neV, 1000*neV,
						15, 15,
						0.01*degree);
      ucn_material[imat]->SetMaterialPropertiesTable(mattbl[imat]);
      materials.push_back(material_name);
      imat++;
      std::cout<<"done."<<std::endl;
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

  G4double worldSizeX = 20.*m;
  G4double worldSizeY = 20.*m;
  G4double worldSizeZ = 20.*m;

  G4Box* solidWorld = new G4Box("World",
                                worldSizeX/2.,worldSizeY/2.,worldSizeZ/2.);

  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,
                                                    fVacuum,
                                                    "World");

  logicWorld->SetUserLimits(g4limit);

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
  char solid_name[100], logical_name[100], physical_name[100];

  G4VisAttributes* VA = new G4VisAttributes(G4Color(0.8,0.8,0.8));
  VA->SetForceSolid(true);

  for (std::map<std::string, std::string>::iterator i = geometryin["GEOMETRY"].begin(); i != geometryin["GEOMETRY"].end(); i++){ // parse STLfile list
    std::istringstream(i->first) >> ID;
    if (ID < 0){
      std::cout << "You defined a solid with ID " << ID << " < 0! IDs have to be larger than zero!\n";
      exit(-1);
    }


    G4Box *btmp = new G4Box("btmp",1*mm,1*mm,1*mm);
    std::istringstream ss(i->second);
    ss >> STLfile >> matname;
    if (ss){
      for (unsigned j = 0; j < materials.size(); j++){
	if (STLfile == "igunored" ||STLfile == "Igunored")break;
	if (matname == "default" ||matname == "Default")break;
	if (matname == materials[j]){
	  if (ID > 1){
	    offset = G4ThreeVector(0, 0, 0);
	    std::cout<<"Reading: "<<STLfile<<std::endl;
	    CADMesh *mesh = new CADMesh((char*)STLfile.c_str(), "STL", m, offset, false);
	    cad_solid[icad] = (G4VSolid*)mesh->TessellatedMesh();
	    sprintf(solid_name,"solid_%d",ID);
	    cad_union[icad] = new G4UnionSolid(solid_name, cad_solid[icad], btmp, 0, G4ThreeVector(5*m, 5*m, 5*m));
	    sprintf(logical_name,"logical_%d",ID);
	    cad_logical[icad] = new G4LogicalVolume(cad_union[icad], ucn_material[j], logical_name, 0, 0, 0);
	    cad_logical[icad]->SetUserLimits(g4limit);
	    cad_logical[icad]->SetVisAttributes(VA);
	    sprintf(physical_name,"physical_%d",ID);
	    cad_physical[icad] = new G4PVPlacement(0, G4ThreeVector(), cad_logical[icad], physical_name, logicWorld, false, 0);
	    icad++;


	  }
	  break;
	}
	else if (j+1 == materials.size()){
	  printf("Material %s used for %s but not defined in geometry.in!\n",matname.c_str(), STLfile.c_str());
	  exit(-1);
	}
      }
    }
    else{
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
    //fField->SetGravityActive(true);

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


void UCNDetectorConstruction::GetMaterial(int imat, std::string name){

  G4NistManager* nistMan = G4NistManager::Instance();

  if(name=="H"||name=="He"||name=="Li"||name=="Be"||name=="B"||
     name=="C"||name=="N"||name=="O"||name=="F"||name=="Ne"||
     name=="Na"||name=="Mg"||name=="Al"||name=="Si"||name=="P"||
     name=="S"||name=="Cl"||name=="Ar"||name=="K"||name=="Ca"||
     name=="Sc"||name=="Ti"||name=="V"||name=="Cr"||name=="Mn"||
     name=="Fe"||name=="Co"||name=="Ni"||name=="Cu"||name=="Zn"||
     name=="Ga"||name=="Ge"||name=="As"||name=="Se"||name=="Br"||
     name=="Kr"||name=="Rb"||name=="Sr"||name=="Y"||name=="Zr"||
     name=="Nb"||name=="Mo"||name=="Tc"||name=="Ru"||name=="Rh"||
     name=="Pd"||name=="Ag"||name=="Cd"||name=="In"||name=="Sn"||
     name=="Sb"||name=="Te"||name=="I"||name=="Xe"||name=="Cs"||
     name=="Ba"||name=="La"||name=="Ce"||name=="Pr"||name=="Nd"||
     name=="Pm"||name=="Sm"||name=="Eu"||name=="Gd"||name=="Tb"||
     name=="Dy"||name=="Ho"||name=="Er"||name=="Tm"||name=="Yb"||
     name=="Lu"||name=="Hf"||name=="Ta"||name=="W"|| name=="Re"||
     name=="Os"||name=="Ir"||name=="Pt"||name=="Au"||name=="Hg"||
     name=="Tl"||name=="Pb"||name=="Bi"||
     //name=="Po"||name=="At"||name=="Rn"||name=="Fr"||
     name=="Ra"||
     //name=="Ac"||
     name=="Th"||name=="Pa"||name=="U"
     //||name=="Np"||name=="Pu"||name=="Am"||
     //name=="Cm"||name=="Bk"||name=="Cf"
     ){
    char g4_name[20]="G4_";
    strcat(g4_name,name.c_str());
    //std::cout<<g4_name<<std::endl;
    ucn_material[imat] = nistMan->FindOrBuildMaterial(g4_name);
    GetElementValues(name);
    return;
  }
  else if(name=="PolishedSteel"){
    ucn_material[imat] = nistMan->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    cohcs=9.4169352; incohcs=1.3042522; scatcs=10.7211872; abscs=3.0528517; fermipot=183.0405; loss=1.04E-04;
    return;
  }
  else if(name=="PE"){
    ucn_material[imat] = nistMan->FindOrBuildMaterial("G4_POLYETHYLENE");
    cohcs=5.005703007; incohcs=11.5357086; scatcs=16.54101522; abscs=0.050797781; fermipot=-8.6553; loss=-2.23E-04;
    return;
  }
  else if(name=="DLC"){
    ucn_material[imat] = new G4Material("DLC", 6, 12.0107*g/mole, 2.8*g/cm3);
    cohcs=5.551; incohcs=0.001; scatcs=5.551; abscs=0.0035; fermipot=243.0964; loss=1.46E-07;
    return;
  }
  else if(name=="DLCLT"){
    ucn_material[imat] = new G4Material("DLCLT", 6, 12.0107*g/mole, 2.8*g/cm3);
    cohcs=5.551; incohcs=0.001; scatcs=5.551; abscs=0.0035; fermipot=243.0964; loss=1.46E-07;
    return;
  }
  else if(name=="UCNdet"){
    ucn_material[imat] = new G4Material("UCNdet", 2.5*g/cm3, 2);
    G4Element* elSi = nistMan->FindOrBuildElement("Si");
    G4Element* elO = nistMan->FindOrBuildElement("O");
    ucn_material[imat] -> AddElement(elSi,1);
    ucn_material[imat] -> AddElement(elO,1);
    cohcs=3.264877149; incohcs=0.002295792; scatcs=3.266746889; abscs=0.080032559; fermipot=90.6365; loss=3.02E-06;
    return;
  }
  else if(name=="NiPLT"){
    ucn_material[imat] = new G4Material("NiPLT", 7.8*g/cm3, 2);
    G4Element* elNi = nistMan->FindOrBuildElement("Ni");
    G4Element* elP = nistMan->FindOrBuildElement("P");
    ucn_material[imat] -> AddElement(elNi,1);
    ucn_material[imat] -> AddElement(elP,1);
    cohcs=11.80105; incohcs=4.42075; scatcs=16.2218; abscs=3.8423; fermipot=212.9579; loss=1.05E-04;
    return;
  }
  else if(name=="NiMo"){
    ucn_material[imat] = new G4Material("NiMo", 9.05*g/cm3, 2);
    G4Element* elNi = nistMan->FindOrBuildElement("Ni");
    G4Element* elMo = nistMan->FindOrBuildElement("Mo");
    ucn_material[imat] -> AddElement(elNi,1);
    ucn_material[imat] -> AddElement(elMo,1);
    cohcs=12.1555; incohcs=4.426; scatcs=16.5815; abscs=4.1885; fermipot=226.7135; loss=1.20E-04;
    return;
  }
  else if(name=="CuBe"){
    ucn_material[imat] = new G4Material("CuBe", 8.36*g/cm3, 2);
    G4Element* elCu = nistMan->FindOrBuildElement("Cu");
    G4Element* elBe = nistMan->FindOrBuildElement("Be");
    ucn_material[imat] -> AddElement(elCu,1);
    ucn_material[imat] -> AddElement(elBe,1);
    cohcs=7.4879; incohcs=0.539036; scatcs=8.022; abscs=3.704552; fermipot=178.8031; loss=1.19E-04;
    return;
  }
  else if(name=="LHe"){
    ucn_material[imat] = new G4Material("LHe", 2, 4.0026*g/mole, 0.1248*g/cm3);
    cohcs=1.34; incohcs=0; scatcs=1.34; abscs=0; fermipot=18.5298; loss=0;//for isopure 4He
    return;
  }
  else if(name=="LHePerf"){
    ucn_material[imat] = new G4Material("LHePerf", 2, 4.0026*g/mole, 0.1248*g/cm3);
    cohcs=1.34; incohcs=0; scatcs=1.34; abscs=0; fermipot=18.5298; loss=0;//for isopure 4He
    return;
  }
  else if(name=="BeO"){
    ucn_material[imat] = new G4Material("BeO", 3.01*g/cm3, 2);
    G4Element* elBe = nistMan->FindOrBuildElement("Be");
    G4Element* elO = nistMan->FindOrBuildElement("O");
    ucn_material[imat] -> AddElement(elBe,1);
    ucn_material[imat] -> AddElement(elO,1);
    cohcs=5.456368569; incohcs=0.00116032; scatcs=5.456368569; abscs=0.002859974; fermipot=256.6665; loss=1.59E-07;
    return;
  }
  else if(name=="SS"){
    //ucn_material[imat] = new G4Material("SS", 8*g/cm3);
    ucn_material[imat] = nistMan->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    cohcs=9.4169352; incohcs=1.3042522; scatcs=10.7211872; abscs=3.0528517; fermipot=183.0405; loss=1.04E-04;
    return;
  }
  else if(name=="dPE"){
    ucn_material[imat] = new G4Material("dPE", 1.071*g/cm3, 2);
    G4Element* elC = nistMan->FindOrBuildElement("C");
    G4Element* elD = new G4Element("Deuterium", "D", 1., 2.014*g/mole);
    ucn_material[imat] -> AddElement(elC,2);
    ucn_material[imat] -> AddElement(elD,4);
    cohcs=5.561297234; incohcs=0.515610557; scatcs=6.075656639; abscs=0.002751316; fermipot=209.4173; loss=6.31E-08;
    return;
  }
  else if(name=="dPS"){
    ucn_material[imat] = new G4Material("dPS", 1.146*g/cm3, 2);
    G4Element* elC = nistMan->FindOrBuildElement("C");
    G4Element* elD = new G4Element("Deuterium", "D", 1., 2.014*g/mole);
    ucn_material[imat] -> AddElement(elC,8);
    ucn_material[imat] -> AddElement(elD,8);
    cohcs=5.55688801; incohcs=0.295256889; scatcs=5.851001289; abscs=0.003071899; fermipot=170.7351; loss=8.39E-08;
    return;
  }
  else{
    std::cout<<"G4Material corresponding to "<<name.c_str()<<" is not defined in UCNDetectorConstruction.cc"<<std::endl;
    exit(-1);
  }
}


void UCNDetectorConstruction::GetElementValues(std::string name){
  if(name=="H"){a=1.00794; cohcs=1.7568; incohcs=80.26; scatcs=82.02; abscs=0.3326; density=0.00008375*g/cm3; fermipot=-0.0487; loss=-2.4719E-05; return;}
  else if(name=="He"){a=4.002602; cohcs=1.34; incohcs=0; scatcs=1.34; abscs=0.00747; density=0.0001663*g/cm3; fermipot=0.0213; loss=6.36747E-07; return;}
  else if(name=="Li"){a=6.941; cohcs=0.454; incohcs=0.92; scatcs=1.37; abscs=70.5; density=0.534*g/cm3; fermipot=-22.9351; loss=-0.010310961; return;}
  else if(name=="Be"){a=9.0121822; cohcs=7.63; incohcs=0.0018; scatcs=7.63; abscs=0.0076; density=1.848*g/cm3; fermipot=250.6328; loss=2.71106E-07; return;}
  else if(name=="B"){a=10.811; cohcs=3.54; incohcs=1.7; scatcs=5.24; abscs=767; density=2.37*g/cm3; fermipot=182.3001; loss=0.040214543; return;}
  else if(name=="C"){a=12.0107; cohcs=5.551; incohcs=0.001; scatcs=5.551; abscs=0.0035; density=1.7*g/cm3; fermipot=147.5942; loss=1.46343E-07; return;}
  else if(name=="N"){a=14.0067; cohcs=11.01; incohcs=0.5; scatcs=11.51; abscs=1.9; density=0.001165*g/cm3; fermipot=0.1222; loss=5.64081E-05; return;}
  else if(name=="O"){a=15.9994; cohcs=4.232; incohcs=0.0008; scatcs=4.232; abscs=0.00019; density=0.001332*g/cm3; fermipot=0.0758; loss=9.09839E-09; return;}
  else if(name=="F"){a=18.99840322; cohcs=4.017; incohcs=0.0008; scatcs=4.018; abscs=0.0096; density=0.00158*g/cm3; fermipot=0.0738; loss=4.71823E-07; return;}
  else if(name=="Ne"){a=20.1797; cohcs=2.62; incohcs=0.008; scatcs=2.628; abscs=0.039; density=0.0008385*g/cm3; fermipot=0.0298; loss=2.37352E-06; return;}
  else if(name=="Na"){a=22.98976928; cohcs=1.66; incohcs=1.62; scatcs=3.28; abscs=0.53; density=0.971*g/cm3; fermipot=24.0558; loss=4.05726E-05; return;}
  else if(name=="Mg"){a=24.305; cohcs=3.631; incohcs=0.08; scatcs=3.71; abscs=0.063; density=1.74*g/cm3; fermipot=60.3755; loss=3.25706E-06; return;}
  else if(name=="Al"){a=26.98153863; cohcs=1.495; incohcs=0.0082; scatcs=1.503; abscs=0.231; density=2.699*g/cm3; fermipot=54.1325; loss=1.86115E-05; return;}
  else if(name=="Si"){a=28.0855; cohcs=2.163; incohcs=0.004; scatcs=2.167; abscs=0.171; density=2.33*g/cm3; fermipot=54.0078; loss=1.14526E-05; return;}
  else if(name=="P"){a=30.97376163; cohcs=3.307; incohcs=0.005; scatcs=3.312; abscs=0.172; density=2.2*g/cm3; fermipot=57.1708; loss=9.31697E-06; return;}
  else if(name=="S"){a=32.065; cohcs=1.0186; incohcs=0.007; scatcs=1.026; abscs=0.53; density=2*g/cm3; fermipot=27.8622; loss=5.17311E-05; return;}
  else if(name=="Cl"){a=35.453; cohcs=11.5257; incohcs=5.3; scatcs=16.8; abscs=33.5; density=0.002995*g/cm3; fermipot=0.1269; loss=0.000972028; return;}
  else if(name=="Ar"){a=39.948; cohcs=0.458; incohcs=0.225; scatcs=0.683; abscs=0.675; density=0.001662*g/cm3; fermipot=0.0125; loss=9.82565E-05; return;}
  else if(name=="K"){a=39.0983; cohcs=1.69; incohcs=0.27; scatcs=1.96; abscs=2.1; density=0.862*g/cm3; fermipot=12.6953; loss=0.000159007; return;}
  else if(name=="Ca"){a=40.078; cohcs=2.78; incohcs=0.05; scatcs=2.83; abscs=0.43; density=1.55*g/cm3; fermipot=28.5202; loss=2.54234E-05; return;}
  else if(name=="Sc"){a=44.9559119; cohcs=19; incohcs=4.5; scatcs=23.5; abscs=27.5; density=2.989*g/cm3; fermipot=128.2093; loss=0.000621791; return;}
  else if(name=="Ti"){a=47.867; cohcs=1.485; incohcs=2.87; scatcs=4.35; abscs=6.09; density=4.54*g/cm3; fermipot=-51.1628; loss=-0.000492238; return;}
  else if(name=="V"){a=50.9415; cohcs=0.0184; incohcs=5.08; scatcs=5.1; abscs=5.08; density=6.11*g/cm3; fermipot=-7.1964; loss=-0.003691556; return;}
  else if(name=="Cr"){a=51.9961; cohcs=1.66; incohcs=1.83; scatcs=3.49; abscs=3.05; density=7.18*g/cm3; fermipot=78.7565; loss=0.000233163; return;}
  else if(name=="Mn"){a=54.9380451; cohcs=1.75; incohcs=0.4; scatcs=2.15; abscs=13.3; density=7.44*g/cm3; fermipot=-79.2569; loss=-0.000990847; return;}
  else if(name=="Fe"){a=55.845; cohcs=11.22; incohcs=0.4; scatcs=11.62; abscs=2.56; density=7.874*g/cm3; fermipot=209.0602; loss=7.52786E-05; return;}
  else if(name=="Co"){a=58.933195; cohcs=0.779; incohcs=4.8; scatcs=5.6; abscs=37.18; density=8.9*g/cm3; fermipot=59.0008; loss=0.004149289; return;}
  else if(name=="Ni"){a=58.6934; cohcs=13.3; incohcs=5.2; scatcs=18.5; abscs=4.49; density=8.902*g/cm3; fermipot=245.1117; loss=0.000121136; return;}
  else if(name=="Cu"){a=63.546; cohcs=7.485; incohcs=0.55; scatcs=8.03; abscs=3.78; density=8.96*g/cm3; fermipot=170.7470; loss=0.000136098; return;}
  else if(name=="Zn"){a=65.38; cohcs=4.054; incohcs=0.077; scatcs=4.131; abscs=1.11; density=7.133*g/cm3; fermipot=97.2309; loss=5.43048E-05; return;}
  else if(name=="Ga"){a=69.723; cohcs=6.675; incohcs=0.16; scatcs=6.83; abscs=2.75; density=5.904*g/cm3; fermipot=96.8294; loss=0.000104855; return;}
  else if(name=="Ge"){a=72.64; cohcs=8.42; incohcs=0.18; scatcs=8.6; abscs=2.2; density=5.323*g/cm3; fermipot=94.1083; loss=7.46909E-05; return;}
  else if(name=="As"){a=74.9215965; cohcs=5.44; incohcs=0.06; scatcs=5.5; abscs=4.5; density=5.73*g/cm3; fermipot=78.9591; loss=0.000190042; return;}
  else if(name=="Se"){a=78.96; cohcs=7.98; incohcs=0.32; scatcs=8.3; abscs=11.7; density=4.5*g/cm3; fermipot=71.2676; loss=0.000407935; return;}
  else if(name=="Br"){a=79.904; cohcs=5.8; incohcs=0.1; scatcs=5.9; abscs=6.9; density=0.007072*g/cm3; fermipot=0.0944; loss=0.000282178; return;}
  else if(name=="Kr"){a=83.798; cohcs=7.67; incohcs=0.01; scatcs=7.68; abscs=25; density=0.003478*g/cm3; fermipot=0.0509; loss=0.000889514; return;}
  else if(name=="Rb"){a=85.4678; cohcs=6.32; incohcs=0.5; scatcs=6.8; abscs=0.38; density=1.532*g/cm3; fermipot=19.9403; loss=1.48936E-05; return;}
  else if(name=="Sr"){a=87.62; cohcs=6.19; incohcs=0.06; scatcs=6.25; abscs=1.28; density=2.54*g/cm3; fermipot=31.9298; loss=5.06683E-05; return;}
  else if(name=="Y"){a=88.9058483; cohcs=7.55; incohcs=0.15; scatcs=7.7; abscs=1.28; density=4.469*g/cm3; fermipot=61.1238; loss=4.58957E-05; return;}
  else if(name=="Zr"){a=91.224; cohcs=6.44; incohcs=0.02; scatcs=6.46; abscs=0.185; density=6.506*g/cm3; fermipot=80.1210; loss=7.17997E-06; return;}
  else if(name=="Nb"){a=92.9063781; cohcs=6.253; incohcs=0.0024; scatcs=6.255; abscs=1.15; density=8.57*g/cm3; fermipot=102.0938; loss=4.53029E-05; return;}
  else if(name=="Mo"){a=95.96; cohcs=5.67; incohcs=0.04; scatcs=5.71; abscs=2.48; density=10.22*g/cm3; fermipot=112.2109; loss=0.000102629; return;}
  else if(name=="Tc"){a=98; cohcs=5.8; incohcs=0.5; scatcs=6.3; abscs=20; density=11.5*g/cm3; fermipot=125.2014; loss=0.000817306; return;}
  else if(name=="Ru"){a=101.07; cohcs=6.21; incohcs=0.4; scatcs=6.6; abscs=2.56; density=12.41*g/cm3; fermipot=135.4358; loss=0.000101192; return;}
  else if(name=="Rh"){a=102.905504; cohcs=4.34; incohcs=0.3; scatcs=4.6; abscs=144.8; density=12.41*g/cm3; fermipot=111.2600; loss=0.006843131; return;}
  else if(name=="Pd"){a=106.42; cohcs=4.39; incohcs=0.093; scatcs=4.48; abscs=6.9; density=12.02*g/cm3; fermipot=104.7363; loss=0.000324433; return;}
  else if(name=="Ag"){a=107.8682; cohcs=4.407; incohcs=0.58; scatcs=4.99; abscs=63.3; density=10.5*g/cm3; fermipot=90.4467; loss=0.002970291; return;}
  else if(name=="Cd"){a=112.411; cohcs=3.04; incohcs=3.46; scatcs=6.5; abscs=2520; density=8.65*g/cm3; fermipot=58.7983; loss=0.143792165; return;}
  else if(name=="In"){a=114.818; cohcs=2.08; incohcs=0.54; scatcs=2.62; abscs=193.8; density=7.31*g/cm3; fermipot=40.6066; loss=0.0132482; return;}
  else if(name=="Sn"){a=118.71; cohcs=4.871; incohcs=0.022; scatcs=4.892; abscs=0.626; density=7.31*g/cm3; fermipot=60.1448; loss=2.79446E-05; return;}
  else if(name=="Sb"){a=121.76; cohcs=3.9; incohcs=0.007; scatcs=3.9; abscs=4.91; density=6.691*g/cm3; fermipot=48.0253; loss=0.000244957; return;}
  else if(name=="Te"){a=127.6; cohcs=4.23; incohcs=0.09; scatcs=4.32; abscs=4.7; density=6.24*g/cm3; fermipot=44.5031; loss=0.000225182; return;}
  else if(name=="I"){a=126.904473; cohcs=3.5; incohcs=0.31; scatcs=3.81; abscs=6.15; density=4.93*g/cm3; fermipot=32.1834; loss=0.000323672; return;}
  else if(name=="Xe"){a=131.293; cohcs=2.96; incohcs=0; scatcs=0; abscs=23.9; density=0.005485*g/cm3; fermipot=0.0322; loss=0.001349884; return;}
  else if(name=="Cs"){a=132.9054519; cohcs=3.69; incohcs=0.21; scatcs=3.9; abscs=29; density=1.873*g/cm3; fermipot=11.9846; loss=0.001486834; return;}
  else if(name=="Ba"){a=137.327; cohcs=3.23; incohcs=0.15; scatcs=3.38; abscs=1.1; density=3.5*g/cm3; fermipot=20.2744; loss=6.02904E-05; return;}
  else if(name=="La"){a=138.90547; cohcs=8.53; incohcs=1.13; scatcs=9.66; abscs=8.97; density=6.154*g/cm3; fermipot=57.2788; loss=0.000302502; return;}
  else if(name=="Ce"){a=140.116; cohcs=2.94; incohcs=0.001; scatcs=2.94; abscs=0.63; density=6.657*g/cm3; fermipot=36.0798; loss=3.61709E-05; return;}
  else if(name=="Pr"){a=140.9076528; cohcs=2.64; incohcs=0.015; scatcs=2.66; abscs=11.5; density=6.71*g/cm3; fermipot=34.2201; loss=0.000697744; return;}
  else if(name=="Nd"){a=144.242; cohcs=7.43; incohcs=9.2; scatcs=16.6; abscs=50.5; density=6.9*g/cm3; fermipot=57.7181; loss=0.001824856; return;}
  else if(name=="Pm"){a=145; cohcs=20; incohcs=1.3; scatcs=21.3; abscs=168.4; density=7.22*g/cm3; fermipot=98.4392; loss=0.003713943; return;}
  else if(name=="Sm"){a=150.36; cohcs=0.422; incohcs=39; scatcs=39; abscs=5922; density=7.46*g/cm3; fermipot=6.2277; loss=2.057036785; return;}
  else if(name=="Eu"){a=151.964; cohcs=6.57; incohcs=2.5; scatcs=9.2; abscs=4530; density=5.243*g/cm3; fermipot=39.0845; loss=0.174351082; return;}
  else if(name=="Gd"){a=157.25; cohcs=29.3; incohcs=151; scatcs=180; abscs=49700; density=7.9*g/cm3; fermipot=51.2364; loss=2.124744306; return;}
  else if(name=="Tb"){a=158.9253468; cohcs=6.84; incohcs=0.004; scatcs=6.84; abscs=23.4; density=8.229*g/cm3; fermipot=59.9568; loss=0.000881096; return;}
  else if(name=="Dy"){a=162.5; cohcs=35.9; incohcs=54.4; scatcs=90.3; abscs=994; density=8.55*g/cm3; fermipot=139.5173; loss=0.016344187; return;}
  else if(name=="Ho"){a=164.9303221; cohcs=8.06; incohcs=0.36; scatcs=8.42; abscs=64.7; density=8.795*g/cm3; fermipot=67.0188; loss=0.002244582; return;}
  else if(name=="Er"){a=167.259; cohcs=7.63; incohcs=1.1; scatcs=8.7; abscs=159; density=9.066*g/cm3; fermipot=66.2510; loss=0.005671831; return;}
  else if(name=="Tm"){a=168.9342133; cohcs=6.28; incohcs=0.1; scatcs=6.38; abscs=100; density=9.321*g/cm3; fermipot=61.2058; loss=0.003930468; return;}
  else if(name=="Yb"){a=173.054; cohcs=19.42; incohcs=4; scatcs=23.4; abscs=34.8; density=6.73*g/cm3; fermipot=75.8461; loss=0.000777986; return;}
  else if(name=="Lu"){a=174.9668; cohcs=6.53; incohcs=0.7; scatcs=7.2; abscs=74; density=9.84*g/cm3; fermipot=63.6214; loss=0.002852069; return;}
  else if(name=="Hf"){a=178.49; cohcs=7.6; incohcs=2.6; scatcs=10.2; abscs=104.1; density=13.31*g/cm3; fermipot=90.0914; loss=0.003756848; return;}
  else if(name=="Ta"){a=180.94788; cohcs=6; incohcs=0.01; scatcs=6.01; abscs=20.6; density=16.65*g/cm3; fermipot=99.7625; loss=0.000828424; return;}
  else if(name=="W"){a=183.84; cohcs=2.97; incohcs=1.63; scatcs=4.6; abscs=18.3; density=19.3*g/cm3; fermipot=80.0538; loss=0.001046354; return;}
  else if(name=="Re"){a=186.207; cohcs=10.6; incohcs=0.9; scatcs=11.5; abscs=89.7; density=21.02*g/cm3; fermipot=162.9495; loss=0.00270937; return;}
  else if(name=="Os"){a=190.23; cohcs=14.4; incohcs=0.3; scatcs=14.7; abscs=16; density=22.57*g/cm3; fermipot=199.1888; loss=0.000415528; return;}
  else if(name=="Ir"){a=192.217; cohcs=14.1; incohcs=0; scatcs=14; abscs=425; density=22.42*g/cm3; fermipot=193.9895; loss=0.011141578; return;}
  else if(name=="Pt"){a=195.084; cohcs=11.58; incohcs=0.13; scatcs=11.71; abscs=10.3; density=21.45*g/cm3; fermipot=165.6172; loss=0.000298146; return;}
  else if(name=="Au"){a=196.9665687; cohcs=7.32; incohcs=0.43; scatcs=7.75; abscs=98.65; density=19.32*g/cm3; fermipot=117.4269; loss=0.003592826; return;}
  else if(name=="Hg"){a=200.59; cohcs=20.24; incohcs=6.6; scatcs=26.8; abscs=372.3; density=13.55*g/cm3; fermipot=134.5205; loss=0.008151295; return;}
  else if(name=="Tl"){a=204.3833; cohcs=9.678; incohcs=0.21; scatcs=9.89; abscs=3.43; density=11.72*g/cm3; fermipot=78.9600; loss=0.000108608; return;}
  else if(name=="Pb"){a=207.2; cohcs=11.115; incohcs=0.003; scatcs=11.118; abscs=0.171; density=11.35*g/cm3; fermipot=80.8339; loss=5.05244E-06; return;}
  else if(name=="Bi"){a=208.9803987; cohcs=9.148; incohcs=0.0084; scatcs=9.156; abscs=0.0338; density=9.747*g/cm3; fermipot=62.4374; loss=1.10085E-06; return;}
  else if(name=="Ra"){a=226; cohcs=13; incohcs=0; scatcs=13; abscs=12.8; density=5*g/cm3; fermipot=34.7128; loss=0.000355692; return;}
  else if(name=="Th"){a=232.03806; cohcs=13.36; incohcs=0; scatcs=13.36; abscs=7.37; density=11.72*g/cm3; fermipot=81.7063; loss=0.000198643; return;}
  else if(name=="Pa"){a=231.035884; cohcs=10.4; incohcs=0.1; scatcs=10.5; abscs=200.6; density=15.37*g/cm3; fermipot=94.9870; loss=0.006125664; return;}
  else if(name=="U"){a=238.02891; cohcs=8.903; incohcs=0.005; scatcs=8.908; abscs=7.57; density=18.95*g/cm3; fermipot=105.1394; loss=0.000249921; return;}
}
