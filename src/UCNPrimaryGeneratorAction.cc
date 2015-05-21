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
UCNPrimaryGeneratorAction::UCNPrimaryGeneratorAction(TConfig GEOMIN, TConfig PARTIN)
{
  geometryin=GEOMIN;
  particlein=PARTIN;
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


  t = G4UniformRand() * ActiveTime * s;
  fParticleGun->SetParticleTime(t);
  RandomPointInSourceVolume();
  //fParticleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, 0.0));
  fParticleGun->SetParticlePosition(G4ThreeVector(X, Y, Z));
  polarization = DicePolarisation(ParticleName);
  fParticleGun->SetParticlePolarization(G4ThreeVector(0,polarization,0));
  particleEnergy = Spectrum(ParticleName) * eV;
  fParticleGun->SetParticleEnergy(particleEnergy);
  AngularDist(ParticleName, phi, theta);
  pz = std::sin(theta)*std::cos(phi);
  px = std::sin(theta)*std::sin(phi);
  py = std::cos(theta);  
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px,py,pz));
  fParticleGun->GeneratePrimaryVertex(anEvent);

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

  std::cout<<"Source type is "<<sourcemode<<std::endl;

  if (sourcemode == "boxvolume"){
    sourceconf >> x_min >> x_max >> y_min >> y_max >> z_min >> z_max >> ActiveTime >> PhaseSpaceWeighting;
  }
  else if (sourcemode == "cylvolume"){
    sourceconf >> r_min >> r_max >> phi_min >> phi_max >> z_min >> z_max >> ActiveTime >> PhaseSpaceWeighting;
    phi_max*=conv;
    phi_min*=conv;
  }
  else if (sourcemode == "cylsurface"){
    sourceconf >> r_min >> r_max >> phi_min >> phi_max >> z_min >> z_max >> ActiveTime >> E_normal;
    phi_max*=conv;
    phi_min*=conv;
    topsurf=(r_max*r_max - r_min*r_min)*(phi_max - phi_min)/2;
    insurf=r_min*(phi_max - phi_min)*(z_max - z_min);
    outsurf=r_max*(phi_max - phi_min)*(z_max - z_min);
    totsurf=topsurf*2+insurf+outsurf;
  }
  else if (sourcemode == "STLvolume"){
    sourceconf >> sourcefile >> ActiveTime >> PhaseSpaceWeighting;
    scene = importer.ReadFile(sourcefile,
			      aiProcess_Triangulate           |
			      aiProcess_JoinIdenticalVertices |
			      aiProcess_CalcTangentSpace);    
    aim = scene->mMeshes[0];
    X_min=1e2*m; X_max=-1e2*m;
    Y_min=1e2*m; Y_max=-1e2*m;
    Z_min=1e2*m; Z_max=-1e2*m;
    G4double X_mesh[3], Y_mesh[3], Z_mesh[3];
    for(unsigned int i=0; i < aim->mNumFaces; i++){
      const aiFace& face = aim->mFaces[i];
      for(int j=0;j<3;j++){
	X_mesh[j] = aim->mVertices[face.mIndices[j]].x*m;
	Y_mesh[j] = aim->mVertices[face.mIndices[j]].y*m;
	Z_mesh[j] = aim->mVertices[face.mIndices[j]].z*m;
      }
      Compare(X_max, X_min, X_mesh);
      Compare(Y_max, Y_min, Y_mesh);
      Compare(Z_max, Z_min, Z_mesh);
    }
  }
  else if (sourcemode == "STLsurface"){
    sourceconf >> sourcefile >> ActiveTime >> E_normal;
    scene = importer.ReadFile(sourcefile,
			      aiProcess_Triangulate           |
			      aiProcess_JoinIdenticalVertices |
			      aiProcess_CalcTangentSpace);    
    aim = scene->mMeshes[0];
    totsurf = 0;
    G4double X_mesh[3], Y_mesh[3], Z_mesh[3];
    for(unsigned int i=0; i < aim->mNumFaces; i++){
      const aiFace& face = aim->mFaces[i];
      for(int j=0;j<3;j++){
	X_mesh[j] = aim->mVertices[face.mIndices[j]].x*m;
	Y_mesh[j] = aim->mVertices[face.mIndices[j]].y*m;
	Z_mesh[j] = aim->mVertices[face.mIndices[j]].z*m;
      }
      totsurf += CalcSurf(X_mesh, Y_mesh, Z_mesh);
    }
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
    X = x_min + G4UniformRand()*(x_max - x_min)*cm;
    Y = y_min + G4UniformRand()*(y_max - y_min)*cm;
    Z = z_min + G4UniformRand()*(z_max - z_min)*cm;
  }
  else if (sourcemode == "cylvolume"){
    r = sqrt(G4UniformRand()*(r_max*r_max - r_min*r_min) + r_min*r_min); // weighting because of the volume element and a r^2 probability outwards
    phi_r = phi_min + G4UniformRand()*(phi_max - phi_min);
    X = r*cos(phi_r)*cm;
    Y = r*sin(phi_r)*cm;
    Z = z_min + G4UniformRand()*(z_max - z_min)*cm;
  }
  else if (sourcemode == "cylsurface"){
    double randsurf=G4UniformRand()*totsurf;

    if(randsurf<outsurf){
      phi_r = phi_min + G4UniformRand()*(phi_max - phi_min);
      X = r_max*cos(phi_r)*cm;
      Y = r_max*sin(phi_r)*cm;
      Z = z_min + G4UniformRand()*(z_max - z_min)*cm;
    }
    else if(randsurf<outsurf+insurf){
      phi_r = phi_min + G4UniformRand()*(phi_max - phi_min);
      X = r_min*cos(phi_r)*cm;
      Y = r_min*sin(phi_r)*cm;
      Z = z_min + G4UniformRand()*(z_max - z_min)*cm;
    }
    else if(randsurf<outsurf+insurf+topsurf){
      r = sqrt(G4UniformRand()*(r_max*r_max - r_min*r_min) + r_min*r_min);
      phi_r = phi_min + G4UniformRand()*(phi_max - phi_min);
      X = r*cos(phi_r)*cm;
      Y = r*sin(phi_r)*cm;
      Z = z_max*cm;
    }
    else{
      r = sqrt(G4UniformRand()*(r_max*r_max - r_min*r_min) + r_min*r_min);
      phi_r = phi_min + G4UniformRand()*(phi_max - phi_min);
      X = r*cos(phi_r)*cm;
      Y = r*sin(phi_r)*cm;
      Z = z_min*cm;
    }
  }
  else if (sourcemode == "STLsurface"){
    G4double randsurf = G4UniformRand()*totsurf;
    G4double tmpsurf = 0;
    G4double X_mesh[3], Y_mesh[3], Z_mesh[3];
    for(unsigned int i=0; i < aim->mNumFaces; i++){
      const aiFace& face = aim->mFaces[i];
      for(int j=0;j<3;j++){
	X_mesh[j] = aim->mVertices[face.mIndices[j]].x*m;
	Y_mesh[j] = aim->mVertices[face.mIndices[j]].y*m;
	Z_mesh[j] = aim->mVertices[face.mIndices[j]].z*m;
      }
      tmpsurf += CalcSurf(X_mesh, Y_mesh, Z_mesh);
      if(randsurf < tmpsurf){
	SetSurfPoint(X_mesh, Y_mesh, Z_mesh);
	break;
      }
    }
    
  }
  else if (sourcemode == "STLvolume"){
    while(1){
      X = X_min+G4UniformRand()*(X_max-X_min);
      Y = Y_min+G4UniformRand()*(Y_max-Y_min);
      Z = Z_min+G4UniformRand()*(Z_max-Z_min);
      if(InSolid(X,Y,Z,Z_min))break;
    }
  }

}

void UCNPrimaryGeneratorAction::Compare(G4double &A_max, G4double &A_min, G4double *A_tmp){
  for(int i=0;i<3;i++){
    if(A_min > A_tmp[i]) A_min = A_tmp[i];
    if(A_max < A_tmp[i]) A_max = A_tmp[i];
  }
}

bool UCNPrimaryGeneratorAction::InSolid(G4double x1, G4double y1, G4double z1, G4double z1_min){
  int ncol = 0;
  G4double X_mesh[3], Y_mesh[3], Z_mesh[3];
  for(unsigned int i=0; i < aim->mNumFaces; i++){
    const aiFace& face = aim->mFaces[i];
    for(int j=0;j<3;j++){
      X_mesh[j] = aim->mVertices[face.mIndices[j]].x*m;
      Y_mesh[j] = aim->mVertices[face.mIndices[j]].y*m;
      Z_mesh[j] = aim->mVertices[face.mIndices[j]].z*m;
    }
    if(Collision(X_mesh, Y_mesh, Z_mesh, x1, y1, z1, z1_min))ncol++;
  }
  if(ncol%2 == 1) return true;
  else            return false;
}


bool UCNPrimaryGeneratorAction::Collision(G4double *X_tmp, G4double *Y_tmp, G4double *Z_tmp, G4double x1, G4double y1, G4double z1, G4double z1_min){
  G4double x2, y2, z2;
  G4double a, b, c, d;// Plane equation: ax+by+cz+d=0
  G4double cp1,cp2,cp3;// Cross products
  x2=x1; y2=y1; z2=z1_min-50*mm;
  a = (Y_tmp[1]-Y_tmp[0])*(Z_tmp[2]-Z_tmp[0])-(Y_tmp[2]-Y_tmp[0])*(Z_tmp[1]-Z_tmp[0]);
  b = (Z_tmp[1]-Z_tmp[0])*(X_tmp[2]-X_tmp[0])-(Z_tmp[2]-Z_tmp[0])*(X_tmp[1]-X_tmp[0]);
  c = (X_tmp[1]-X_tmp[0])*(Y_tmp[2]-Y_tmp[0])-(X_tmp[2]-X_tmp[0])*(Y_tmp[1]-Y_tmp[0]);
  d = -(a*X_tmp[0]+b*Y_tmp[0]+c*Z_tmp[0]);
  if((a*x1+b*y1+c*z1+d)*(a*x2+b*y2+c*z2+d)>=0)
    return false;//Check if the two points are in difference sides of the triangle plane
  cp1=(x1-X_tmp[0])*(Y_tmp[1]-Y_tmp[0])-(y1-Y_tmp[0])*(X_tmp[1]-X_tmp[0]);
  cp2=(x1-X_tmp[1])*(Y_tmp[2]-Y_tmp[1])-(y1-Y_tmp[1])*(X_tmp[2]-X_tmp[1]);
  cp3=(x1-X_tmp[2])*(Y_tmp[0]-Y_tmp[2])-(y1-Y_tmp[2])*(X_tmp[0]-X_tmp[2]);
  if(!(cp1<0&&cp2<0&&cp3<0)&&!(cp1>0&&cp2>0&&cp3>0))
    return false;//Check if the interpolated point between the two points is in the triangle.    
  else
    return true;
}

G4double UCNPrimaryGeneratorAction::CalcSurf(G4double *X_tmp, G4double *Y_tmp, G4double *Z_tmp){
  G4double xg, yg, zg, surf;
  xg = (Y_tmp[1]-Y_tmp[0])*(Z_tmp[2]-Z_tmp[0])-(Z_tmp[1]-Z_tmp[0])*(Y_tmp[2]-Y_tmp[0]);
  yg = (Z_tmp[1]-Z_tmp[0])*(X_tmp[2]-X_tmp[0])-(X_tmp[1]-X_tmp[0])*(Z_tmp[2]-Z_tmp[0]);
  zg = (X_tmp[1]-X_tmp[0])*(Y_tmp[2]-Y_tmp[0])-(Y_tmp[1]-Y_tmp[0])*(X_tmp[2]-X_tmp[0]);
  surf = sqrt(xg*xg + yg*yg + zg*zg)/2;
  return surf;
}

void UCNPrimaryGeneratorAction::SetSurfPoint(G4double *X_tmp, G4double *Y_tmp, G4double *Z_tmp){
  G4double t1=sqrt(G4UniformRand());
  G4double t2=1-t1;
  G4double s1=G4UniformRand();
  G4double s2=1-s1;
  G4double xtmp[2],ytmp[2],ztmp[2];
  xtmp[0]=X_tmp[0]*t1+X_tmp[2]*t2;
  ytmp[0]=Y_tmp[0]*t1+Y_tmp[2]*t2;
  ztmp[0]=Z_tmp[0]*t1+Z_tmp[2]*t2;
  xtmp[1]=X_tmp[1]*t1+X_tmp[2]*t2;
  ytmp[1]=Y_tmp[1]*t1+Y_tmp[2]*t2;
  ztmp[1]=Z_tmp[1]*t1+Z_tmp[2]*t2;
  X=xtmp[0]*s1+xtmp[1]*s2;
  Y=ytmp[0]*s1+ytmp[1]*s2;
  Z=ztmp[0]*s1+ztmp[1]*s2;
}
