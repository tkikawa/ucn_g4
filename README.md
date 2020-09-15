UCN simulation code with GEANT4
========

This is the Monte Calro simulation code of UCN (ultra-cold neutron) with GEANT4 developed for [the TRIUMF UCN experiment](http://www.triumf.ca/ucn). It has functions equivalent to [PENTrack](https://github.com/wschreyer/PENTrack) and intput and output format is the same as PENTrack.


## 1. External libraries

### GEANT4

[GEANT4](http://geant4.cern.ch/) is a toolkit for the simulation of the passage of particles through matter developed by CERN. GEANT4.10.01.p01 has been tested. If the GEANT4 library is used, the following reference publication should be quoted:  
[S. Agostinelli et al. (GEANT4 Collaboration), Nucl. Instrum. Meth. B 506, 250 (2003)](http://www.sciencedirect.com/science/article/pii/S0168900203013688).

### ASSIMP

[ASSIMP](http://assimp.sourceforge.net/) is a portable open source library to import [various well-known 3D model formats](http://assimp.sourceforge.net/main_features_formats.html) in a uniform manner. ASSIMP v.3.1.1
has been tested.

### muparser

[muparser](http://muparser.beltoforion.de/) is a fast formula parser and is used to interpret energy distributions etc. given by the user in particle.in. muparser v.2.2.4 has been tested.


### ALGLIB

[ALGLIB](http://www.alglib.net) is used to do 1D and 2D interpolation for field calculations. It is included in the repository.


### libtricubic

[Lekien and Marsden](http://dx.doi.org/10.1002/nme.1296) developed a tricubic interpolation method in three dimensions. It is included in the repository.


## 2. Installation

For the installtion, [cmake](http://www.cmake.org/) is used. The cmake version newer than 2.8 must be used. (cmake v.2.6 is normally installed in Scientific Linux 6.)


### Installation of GEANT4

- Download the codes from [GEANT4 website](http://geant4.cern.ch/support/download.shtml) and decompress it.
- export G4INSTALL=/usr/local/geant4/ (directory where GEANT4 will be installed.)  
- mkdir build (Any name of the directory is OK.)  
- cd build (You must run cmake in other directory where CMakeLists.txt exists.)  
- cmake ../ -DCMAKE_INSTALL_PREFIX=/usr/local/geant4/ (directory where GEANT4 will be installed.)  
  -DXERCESC_ROOT_DIR=/usr/include (directory where xercec-c is installed.)  
  -DGEANT4_INSTALL_DATA=ON (install the additional data.)  
  -DGEANT4_USE_OPENGL_X11=ON (activate OpenGL.)  
  $G4INSTALL  
- make (It takes a lot of time. If you'd like to compile with more than one core, input "make -j {number of cores}".)  
- make install

### Installation of ASSIMP

- Download the codes from [ASSIMP website](http://assimp.sourceforge.net/main_downloads.html) and decompress it.  
- mkdir build  
- cd build  
- cmake ../ -DCMAKE_INSTALL_PREFIX=/usr/local/assimp/ (directory where ASSIMP will be installed.)  
- make  
- make install  

### Installation of muparser

- Download the codes from [muparser website](http://muparser.beltoforion.de/mup_download.html) and decompress it.  
- ./configure --prefix=/usr/local/muparser (directory where muparser will be installed.)  
- make  
- make install  

### Installation of this simulation code
- (Edit the directory settings in CMakeLists.txt if you installed the abeve libraries in other directories than defaults.)  
- mkdir build  
- cd build  
- cmake ../  
- make  


## 3. Run the simulation

### Usage
- UCNsim {jobnumber} {path/to/in/files} {path/to/out/files} [options]
- Default settings: jobnumber=0, path/to/in/files="./in", path/to/out/files="./out"

### Options
- -m : set the macro
- -s : use the jobnumber as the seed
- -t : set the number of threads

## 4. Defining your experiment

### Geometry

Geometry can be imported using [binary STL files] (http://en.wikipedia.org/wiki/STL_%28file_format%29).
STL files use a list of triangles to describe 3D-surfaces. They can be created with most CAD software, e.g. Solidworks via "File - Save As..." using the STL file type. Solidworks makes sure that the surface normals always point outside a solid body. This is used to check if points lie inside a solid object.
Note the "Options..." button in the Save dialog, there you can change format (only binary supported), resolution, coordinate system etc.

Gravity acts in negative z-direction, so choose your coordinate system accordingly.

Do not choose a too high resolution. Usually, setting the deviation tolerance to a value small enough to represent every detail in your geometry and leaving angle tolerance at the most coarse setting is good enough. Low angle tolerance quickly increases triangle count and degeneracy.

If you want to export several parts of a Solidworks assembly you can do the following:

1. Select the part(s) to be exported and rightclick.
2. Select "Invert selection" and then use "Suppress" on all other parts.
3. Now you can save that single part as STL (make sure the option "Do not translate STL output data to positive space" is checked and to use the same coordinate system for every part or else the different parts will not fit together).
4. You can check the positioning of the parts with e.g. [MeshLab](http://meshlab.sourceforge.net/), SolidView, Minimagics, Solidworks...

### Fields

Magnetic and electric fields (rotationally symmetric 2D and 3D) can be included from text-based field maps
(right now, only "[Vectorfields OPERA](https://www.cobham.com/about-cobham/aerospace-and-security/about-us/antenna-systems/specialist-technical-services-and-software/products-and-services/design-simulation-software/opera.aspx)" maps in cgs units are supported).
You can also define analytic fields from straight, finite conductors.

### Particle sources

Particle sources can be defined using STL files or manual parameter ranges. Particle spectra and velocity distributions can also be conviniently defined in the particle.in file.



## 5. Output

Output files are separated by particle type, (e.g. electron, neutron and proton) and type of output (endlog, tracklog, ...). Output files are only created if particles of the specific type are simulated and can also be completely disabled for each particle type individually by adding corresponding variables in 'particle.in'. All output files are tables with space separated columns; the first line contains the column name.

Types of output: endlog, tracklog, hitlog, snapshotlog, spinlog.

### Endlog

The endlog keeps track of the starting and end parameters of the particles simulated. In the endlog, you get the following parameters:

- jobnumber: corresponds to job number of the run (command line parameter). Only used when running multiple instances simultaneously.
- particle: number which corresponds to particle being simulated 
- tstart: time at which particle is created in the simulation [s]
- xstart, ystart, zstart: coordinate from which the particle starts [m]
- vxstart, vystart,vzstart: velocity with which the particle starts [m/s]
- polstart: initial polarization of the particle (-1,1)
- Hstart: initial total energy of particle [eV]
- Estart: initial kinetic energy of the particle [eV]
- Bstart: magnetic field at starting point [T]
- Ustart: electric potential at starting point [V]
- solidstart: number of geometry part the particle started in, see geometry.in
- tend: time at which particle simulation is stopped [s]
- xend, yend, zend: coordinate at which the particle ends [m]
- vxend, vyend, vzend: velocity with which the particle ends [m/s]
- polend: final polarization of the particle (-1,1)
- Hend: final total energy of particle [eV]
- Eend: final kinetic energy of the particle [eV]
- Bend: magnetic field at stopping point [T]
- Uend: electric potential at stopping point [V]
- solidend: number of geometry part the particle stopped in, defined in geometry.in
- stopID: code which identifies why the particle was stopped
  - 0: not categorized
  - -1: did not finish (reached max. simulation time)
  - -2: hit outer boundaries
  - -4: decayed
  - 1: absorbed in bulk material (see solidend)
  - 2: absorbed on total reflection on surface (see solidend)
- NSpinflip: number of spin flips that the particle underwent during simulation
- spinflipprob: probability that the particle has undergone a spinflip, calculated by bruteforce integration of the Bloch equation
- ComputingTime: Time which is required to compute and track the particle [s]
- Nhit: number of times particle hit a geometry surface (e.g. wall of guidetube)
- Nstep: number of steps that it took to simulate particle
- trajlength: the total length of the particle trajectory from creation to finish [m]
- Hmax: the maximum total energy that the particle had during trajectory [eV]

### Snapshotlog

Switching on snapshotlog in "particle.in" will output the particle parameters at additional 'snapshot times' (also defined in "particle.in") in the snapshotlog. It contains the same data fields as the endlog.

### Tracklog

The tracklog contains the complete trajectory of the particle, with following parameters:

- jobnumber: corresponds to job number of the run (command line parameter). Only used when running multiple instances simultaneously.
- particle: number which corresponds to particle being simulated 
- polarization: the polarization of the particle (-1,1)
- t: time [s]
- x,y,z coordinates [m]
- vx,vy,vz: velocity [m/s]
- H: total energy [eV]
- E: kinetic energy [eV]
- Bx: x-component of magnetic field at coordinates [T]
- dBxdx, dBxdy, dBxdz: gradient of x-component of magnetic field at coordinates [T/m]
- By: y-component of magnetic field at coordinates[T]
- dBydx, dBydy, dBydz: gradient of y-component of magnetic field at coordinates [T/m]
- Bz: z-component of magnetic field at coordinates[T]
- dBzdx, dBzdy, dBzdz: gradient of z-component of magnetic field at coordinates [T/m]
- Ex, Ey, Ez: X, Y, and Z component of electric field at coordinates [V/m]
- V: electric potential at coordinates [V]

### Hitlog

This log contains all the points at which particle hits a surface of the experiment geometry. This includes both reflections and transmissions.

In the hitlog, you get the following parameters:

- jobnumber: number which corresponds to job number of the run (command line parameter). Only used when running multiple instances simultaneously.
- particle: number which corresponds to particle being simulated
- t: time at which particle hit the surface [s]
- x, y, z: coordinate of the hit [m]
- v1x, v1y,v1z: velocity of the particle before it interacts with the new geometry/surface [m/s]
- pol1: polarization of particle before it interacts with the new geometry/surface (-1,1)
- v2x, v2y,v2z: velocity of the particle after it interacted with the new geometry/surface [m/s]
- pol2: polarization of particle after it interacted with the new geometry/surface (-1,1)
- nx,ny,nz: normal to the surface of the geometry that the particle hits
- solid1: ID number of the geometry that the particle starts in
- solid2: ID number of the geometry that the particle hits

### Spinlog

If bruteforce integration of particle spin precession is active, it will be logged into the spinlog. In the spinlog, you get the following parameters:

- t: time [s]
- Babs: absolute magnetic field [T]
- Polar: projection of Bloch spin vector onto magnetic field direction vector (Bloch vector has dimensionless length 0.5, so this quantity is also dimensionless)
- logPolar: logarithm of above value [dimensionless]
- Ix, Iy, Iz: x, y and z components of the Bloch vector (times 2) [dimensionless]
- Bx, By, Bz: x, y and z components of magnetic field direction vector (B[i]/Babs) [dimensionless]

### Writing Output files to ROOT readable files

merge_all.c: A [ROOT](http://root.cern.ch) script that writes all out-files (or those, whose filename matches some pattern) into a single ROOT file containing trees for each log- and particle-type.
