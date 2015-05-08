#include "UCNPrimaryGeneratorAction.hh"

#include "UCNDetectorConstruction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//UCNPrimaryGeneratorAction::UCNPrimaryGeneratorAction(void)
UCNPrimaryGeneratorAction::UCNPrimaryGeneratorAction(int simcount, int secondaries, TConfig GEOMIN, TConfig PARTIN)
{
  geometryin=GEOMIN;
  particlein=PARTIN;
  ncount=simcount;
  ParticleName="neutron";
  ConfigInit(particlein);
  TSource(geometryin);

  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle(ParticleName);
  fParticleGun->SetParticleDefinition(particle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UCNPrimaryGeneratorAction::~UCNPrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UCNPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event


  for(int imc=0; imc<ncount; imc++){

    t = G4UniformRand() * ActiveTime;
    fParticleGun->SetParticleTime(t);
    RandomPointInSourceVolume();
    //fParticleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, 0.0));
    fParticleGun->SetParticlePosition(G4ThreeVector(X*cm, Y*cm, Z*cm));
    polarization = DicePolarisation(ParticleName);
    fParticleGun->SetParticlePolarization(G4ThreeVector(0,polarization,0));
    
    particleEnergy = Spectrum(ParticleName);
    fParticleGun->SetParticleEnergy(particleEnergy);
    AngularDist(ParticleName, phi, theta);
    if (phi > pi/2 && phi < pi) phi = pi-phi;
    
    pz = std::sin(phi)*std::cos(theta);
    px = std::sin(phi)*std::sin(theta);
    py = std::cos(phi);
    
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px,py,pz));
    fParticleGun->GeneratePrimaryVertex(anEvent);

  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
void UCNPrimaryGeneratorAction::ConfigInit(std::map<std::string, std::string> conf){
 
  //float nextsnapshot = -1;
  snapshotlog = false;
  std::istringstream(conf["snapshotlog"]) >> snapshotlog;
  //std::istringstream snapshots(conf["snapshots"]);
  if (snapshotlog){
    //do{
    //snapshots >> nextsnapshot;
    //}while (snapshots.good() && nextsnapshot < x); // find first snapshot time
    std::cout<<"Snapshot is not available.\n";
  }

  tracklog = false;
  std::istringstream(conf["tracklog"]) >> tracklog;
  trackloginterval = 1e-3;
  std::istringstream(conf["trackloginterval"]) >> trackloginterval;

  hitlog = false;
  std::istringstream(conf["hitlog"]) >> hitlog;

  std::istringstream(conf["flipspin"]) >> flipspin;


}
*/

void UCNPrimaryGeneratorAction::ConfigInit(TConfig conf){
  for (TConfig::iterator i = conf.begin(); i != conf.end(); i++){
    TParticleConfig *pconf = &pconfigs[i->first];
    std::istringstream(i->second["Emax"]) >> pconf->Emax;
    std::istringstream(i->second["Emin"]) >> pconf->Emin;
    std::istringstream(i->second["lmax"]) >> pconf->lmax;
    std::istringstream(i->second["polarization"]) >> pconf->polarization;
    std::istringstream(i->second["tau"]) >> pconf->tau;
    std::istringstream(i->second["tmax"]) >> pconf->tmax;
    std::istringstream(i->second["phi_v_min"]) >> pconf->phi_v_min;
    std::istringstream(i->second["phi_v_max"]) >> pconf->phi_v_max;
    std::istringstream(i->second["theta_v_min"]) >> pconf->theta_v_min;
    std::istringstream(i->second["theta_v_max"]) >> pconf->theta_v_max;
    try{
      pconf->spectrum.DefineVar("x", &xvar);
      pconf->spectrum.SetExpr(i->second["spectrum"]);
      pconf->spectrum.DefineFun("ProtonBetaSpectrum", &ProtonBetaSpectrum);
      pconf->spectrum.DefineFun("ElectronBetaSpectrum", &ElectronBetaSpectrum);
      pconf->phi_v.DefineVar("x", &xvar);
      pconf->phi_v.SetExpr(i->second["phi_v"]);
      pconf->theta_v.DefineVar("x", &xvar);
      pconf->theta_v.SetExpr(i->second["theta_v"]);
    }
    catch (mu::Parser::exception_type &exc){
      std::cout << exc.GetMsg();
      exit(-1);
    }
  }
}



void UCNPrimaryGeneratorAction::TSource(TConfig geometryconf){

  sourcemode = geometryconf["SOURCE"].begin()->first; // only first source in geometry.in is read in
  std::istringstream sourceconf(geometryconf["SOURCE"].begin()->second);
  sourceconf >> ParticleName;

  if (sourcemode == "boxvolume"){
    sourceconf >> x_min >> x_max >> y_min >> y_max >> z_min >> z_max >> ActiveTime >> PhaseSpaceWeighting;
  }
  else if (sourcemode == "cylvolume"){
    sourceconf >> r_min >> r_max >> phi_min >> phi_max >> z_min >> z_max >> ActiveTime >> PhaseSpaceWeighting;
  }
  else if (sourcemode == "STLvolume"){
    sourceconf >> sourcefile >> ActiveTime >> PhaseSpaceWeighting;
  }
  else if (sourcemode == "cylsurface"){
    sourceconf >> r_min >> r_max >> phi_min >> phi_max >> z_min >> z_max >> ActiveTime >> E_normal;
  }
  else if (sourcemode == "STLsurface"){
    sourceconf >> sourcefile >> ActiveTime >> E_normal;
  }
  else{
    std::cout << "\nCould not load source """ << sourcemode << """! Did you enter invalid parameters?\n";
    exit(-1);
  }


}




double UCNPrimaryGeneratorAction::Spectrum(const std::string &particlename){
  double y;
  TParticleConfig *pconfig = &pconfigs[particlename];
  for (;;){
    xvar = pconfig->Emin+(pconfig->Emax-pconfig->Emin)*G4UniformRand();
    try{
      y = pconfig->spectrum.Eval();
    }
    catch(mu::Parser::exception_type &exc){
      std::cout << exc.GetMsg();
      exit(-1);
    }
    if (G4UniformRand() < y)
      return xvar;
  }
  return 0;
}





void UCNPrimaryGeneratorAction::AngularDist(const std::string &particlename, double &phi_v, double &theta_v){
  double y;
  TParticleConfig *pconfig = &pconfigs[particlename];
  for (;;){
    xvar = pconfig->phi_v_min+(pconfig->phi_v_max-pconfig->phi_v_min)*G4UniformRand();
    try{
      y = pconfig->phi_v.Eval();
    }
    catch(mu::Parser::exception_type &exc){
      std::cout << exc.GetMsg();
      exit(-1);
    }
    if (G4UniformRand() < y){
      phi_v = xvar;
      break;
    }
  }
  for (;;){
    xvar = pconfig->theta_v_min+(pconfig->theta_v_max-pconfig->theta_v_min)*G4UniformRand();
    try{
      y = pconfig->theta_v.Eval();
    }
    catch(mu::Parser::exception_type &exc){
      std::cout << exc.GetMsg();
      exit(-1);
    }
    if (G4UniformRand() < y){
      theta_v = xvar;
      break;
    }
  }
}

int UCNPrimaryGeneratorAction::DicePolarisation(const std::string &particlename){
  int p = pconfigs[particlename].polarization;
  if (p == 0){
    if(G4UniformRand() < 0.5)
      return -1;
    else
      return 1;
  }
  else
    return p;
}

void UCNPrimaryGeneratorAction::RandomPointInSourceVolume(){

  double r, phi_r;

  if (sourcemode == "boxvolume"){
    X = x_min + G4UniformRand()*(x_max - x_min);
    Y = y_min + G4UniformRand()*(y_max - y_min);
    Z = z_min + G4UniformRand()*(z_max - z_min);
  }
  else if (sourcemode == "cylvolume"){
    r = sqrt(G4UniformRand()*(r_max*r_max - r_min*r_min) + r_min*r_min); // weighting because of the volume element and a r^2 probability outwards
    phi_r = phi_min + G4UniformRand()*(phi_max - phi_min);
    X = r*cos(phi_r);
    Y = r*sin(phi_r);
    Z = z_min + G4UniformRand()*(z_max - z_min);
  }
  else{
    std::cout << "STLsurface and STLvolume are not implemented yet.\n";
    exit(-1);
  }

}
