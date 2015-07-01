#ifndef GLOBALS_H_
#define GLOBALS_H_

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <string>
#include <map>

#define ID_UNKNOWN 0 ///< standard flag for particles
#define ID_NOT_FINISH -1 ///< flag for particles which reached ::StorageTime
#define ID_HIT_BOUNDARIES -2 ///< flag for particles which left bounding box of TParticle::geom
//#define ID_ODEINT_ERROR -3 ///< flag for particles which produced a numerical error during ODe integration
#define ID_DECAYED -4 ///< flag for particles which reached TParticle::tau
//#define ID_INITIAL_NOT_FOUND -5 ///< flag for particles which had a too low total energy to find a initial spot in the source volume
//#define ID_CGAL_ERROR -6 ///< flag for particles which produced an error during geometry collision checks
//#define ID_GEOMETRY_ERROR -7 ///< flag for particles which produced an error while tracking material boundaries along the trajectory
#define ID_ABSORBED_IN_MATERIAL 1 ///< flag for particles that were absorbed inside a material
#define ID_ABSORBED_ON_SURFACE 2 ///< flag for particles that were absorbed on a material surface

#define PARTICLE 1 ///< set particletype in configuration to this value to simulate particles
#define BF_ONLY 3 ///< set particletype in configuration to this value to print out a ramp heating analysis
#define BF_CUT 4 ///< set particletype in configuration to this value to print out a planar slice through electric/magnetic fields
#define GEOMETRY 7 ///< set particletype in configuration to this value to print out a sampling of the geometry

// physical constants
static const long double ele_e = 1.602176487E-19L; ///< elementary charge [C]
static const long double gravconst = 9.80665L; ///< g [m/s]
static const long double conv = 3.1415926535897932384626L/180.L; ///< deg to rad conversion factor
static const long double Mu0 = 4*pi*1e-7L; ///< magnetic permeability [Vs/Am]
static const long double m_n = 1.674927211E-27L/ele_e; ///< neutron mass [eV/c^2]
static const long double m_p = 1.672621637E-27L/ele_e; ///< proton mass [eV/c^2]
static const long double m_e = 9.10938215e-31L/ele_e; ///< electron mass [eV/c^2]
static const long double c_0 = 299792458.L; ///< light speed [m/s]
static const long double hbar = 1.05457266e-34L; ///< planck constant [Js]
static const long double mu_nSI = -0.96623641e-26L;	///< Neutron magnetic moment [J/T]
static const long double gamma_n = -1.83247185e8L; ///< 2*::mu_nSI/::hquer gyromagnetic ratio of neutron [1/Ts]

static const long double lengthconv = 0.01; ///< length conversion factor cgs -> SI [cm -> m]
static const long double Bconv = 1e-4; ///< magnetic field conversion factor cgs -> SI [G -> T]
static const long double Econv = 1e2; ///< electric field conversion factor cgs -> SI [V/cm -> V/m]

// gravity
static const long double gx = 0;
static const long double gy = 0;
//static const long double gz = -gravconst;
static const long double gz = 0;


/**
 * Print progress bar.
 *
 * Prints a point every 2% and a number every 10%
 *
 * @param percentage Progress of current action (0..1)
 * @param lastprint Put in and return last printed percentage
 */
void PrintPercent(double percentage, int &lastprint);


/**
 * Energy distribution of protons from free neutron beta decay (0 < E < 750 eV)
 *
 * From diploma thesis M. Simson.
 *
 * @param E energy
 *
 * @return Returns a probability between 0 and 1 that a proton with energy E is created.
 */
double ProtonBetaSpectrum(double E);


/**
 * Energy distribution of electrons from free neutron decay (0 < E < 782 keV)
 *
 * From "http://hyperphysics.phy-astr.gsu.edu/Hbase/nuclear/beta2.html"
 *
 * @param E energy
 *
 * @return Returns probability between 0 and 1 that an electron with energy E is created.
 */
double ElectronBetaSpectrum(double E);


typedef std::map<std::string, std::map<std::string, std::string> > TConfig; ///< map of sections containing a map of key-value pairs

/**
 * Read variables from *.in file into map.
 *
 * in file is assumed to have structure
 * [section]
 * key arbitrary string
 * key value
 * [section]
 * ...
 *
 * @param inpath Path to in file.
 * @param vars Return TConfig map
 */
void ReadInFile(const char *inpath, TConfig &vars);



class TField{
public:
  /**
   * Add magnetic field at a given position and time.
   *
   * Add field components to the given field matrix B:
   *Bx,dBxdx,dBxdy,dBxdz;
   *By,dBydx,dBydy,dBydz;
   *Bz,dBzdx,dBzdy,dBzdz;
   *Babs,dBdx,dBdy,dBdz;
   * Has to be implemented by all derived field calculation classes.
   *
   * @param x Cartesian x coordinate
   * @param y Cartesian y coordinate
   * @param z Cartesian z coordinate
   * @param t Time
   * @param B Magnetic field component matrix to which the values are added
   */
  virtual void BField (double x, double y, double z, double t, double B[4][4]) = 0;

  /**
   * Add electric field and potential at a given position.
   *
   * Has to be implemented by all derived field calculation classes.
   *
   * @param x Cartesian x coordinate
   * @param y Cartesian y coordinate
   * @param z Cartesian z coordinate
   * @param t Time
   * @param V Return electric potential (!=0 only if a map with potential was loaded)
   * @param Ei Returns electric field vector
   */
  virtual void EField (double x, double y, double z, double t, double &V, double Ei[3]) = 0;

  /**
   * Virtual destructor
   */
  virtual ~TField(){ };
};


#endif /*GLOBALS_H_*/
