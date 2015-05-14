/********************************************/
/*                                          */
/*  GEANT4 simulation code for TRIUMF UCN   */
/*                                          */
/********************************************/

#include <cstdlib>
#include <cstdio>
#include <csignal>
#include <cmath>
#include <ctime>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sys/time.h>


#include "UCNPhysicsList.hh"
#include "UCNDetectorConstruction.hh"
#include "UCNActionInitialization.hh"
#include "UCNGlobals.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "G4UImanager.hh"
#include "G4UIcommand.hh"

#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {
  void PrintUsage() {
    G4cerr << "Usage: " << G4endl;
    G4cerr << " UCNsim <jobnumber> <path/to/in/files> <path/to/out/files> [options]"<< G4endl;
    G4cerr << "Options:"<< G4endl;
    G4cerr << " -m : set the macro"<< G4endl;
    G4cerr << " -s : use the jobnumber as the seed"<< G4endl;
    G4cerr << " -t : set the number of threads (available only for multi-threaded mode)"<< G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Evaluate arguments
  //
  if ( argc > 9 ) {
    PrintUsage();
    return 1;
  }

  G4String macro;
#ifdef G4MULTITHREADED
  G4int nThreads = 0;
#endif

  std::string outpath = "./out";
  std::string inpath  = "./in";
  int jobnumber = 0;
  bool jobseed = false;

  int c = -1;
  int op = 0;
  while ((c = getopt(argc, argv, "m:st:")) != -1) {
    switch(c){
    case 'm':
      macro = optarg;
      break;
    case 's':
      jobseed = true;
      break;
    case 't':
#ifdef G4MULTITHREADED
      nThreads = atoi(optarg);
#endif
      break;
    default:
      PrintUsage();
      return 1;
    }
  }

  while (optind < argc){
    switch(op){
    case 0:
      jobnumber = atoi(argv[optind++]);
      break;
    case 1:
      inpath = argv[optind++];
      break;
    case 2:
      outpath = argv[optind++];
      break;
    default:
      G4cerr << "Too many arguments." << G4endl;
      PrintUsage();
      return 1;
    }
    op++;
  }

  TConfig configin;
  ReadInFile(std::string(inpath + "/config.in").c_str(), configin);
  TConfig geometryin;
  ReadInFile(std::string(inpath + "/geometry.in").c_str(), geometryin);
  TConfig particlein;
  ReadInFile(std::string(inpath + "/particle.in").c_str(), particlein);
  for (TConfig::iterator i = particlein.begin(); i != particlein.end(); i++){
    if (i->first != "all"){
      i->second = particlein["all"];
    }
  }
  ReadInFile(std::string(inpath + "/particle.in").c_str(), particlein);


  int simtype = PARTICLE;
  //int neutdist = 0;
  int simcount = 1;
  double SimTime = 1000;
  int secondary = 1;
  std::istringstream(configin["global"]["simtype"])>> simtype;
  //std::istringstream(configin["global"]["neutdist"])>> neutdist;
  std::istringstream(configin["global"]["simcount"])>> simcount;
  std::istringstream(configin["global"]["simtime"])>> SimTime;
  std::istringstream(configin["global"]["secondaries"])>> secondary;
  simtype=1;
  if(simtype!=1){
    std::cout<<"Only the particle simualtion is available.\n";
    std::cout<<"Set simtype in config.in 1.\n";
    exit(-1);
  }


  // Choose the Random engine
  //
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // Construct the default run manager
  //
#ifdef G4MULTITHREADED 
  G4MTRunManager * runManager = new G4MTRunManager;
  if ( nThreads > 0 ) runManager->SetNumberOfThreads(nThreads);
#else
  G4RunManager * runManager = new G4RunManager;
#endif

  // Seed the random number generator manually
  if(jobseed) G4Random::setTheSeed(jobnumber);

  // Set mandatory initialization classes
  //
  // Detector construction
  runManager->SetUserInitialization(new UCNDetectorConstruction(SimTime, geometryin));
  G4cout << "Detector init. OK!" << G4endl;
  // Physics list
  runManager->SetUserInitialization(new UCNPhysicsList());
  G4cout << "PhysicsList init. OK!" << G4endl;
  // User action initialization
  runManager->SetUserInitialization(new UCNActionInitialization(jobnumber, outpath, secondary, geometryin, particlein));
  G4cout << "Action init. OK!" << G4endl;

  // Initialize G4 kernel
  runManager->Initialize();

#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
  G4cout <<"visualization init. OK!" << G4endl;
#endif

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  
  if ( macro.size() ) {
    // batch mode
    G4UIsession * session = 0;
    G4String command = "/control/execute ";
    session = new G4UIterminal(new G4UItcsh);
    UImanager->ApplyCommand(command+macro); 
    session->SessionStart();
    delete session;
  }
  else
  {
    // normal mode
    //runManager->BeamOn(1);
    runManager->BeamOn(simcount);
  }
  
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  G4cout << "All finished. \\(^o^)/" << G4endl;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
