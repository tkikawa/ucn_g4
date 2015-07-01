#ifndef UCNField_H
#define UCNField_H 1

#include "UCNGlobals.hh"
#include "UCN2DField.hh"
#include "UCN3DField.hh"
#include "UCNConductorField.hh"
#include "UCNField.hh"

#include "G4Field.hh"
#include "G4ElectroMagneticField.hh"

//class UCNField : public G4Field
class UCNField : public G4ElectroMagneticField
{
public:
  
  UCNField(std::vector<TField*> FIELDS);
  
  void GetEMFieldValue(const double t, const double x, const double y, const double z, double B[4][4], double *Ei, double &V);
  
  virtual ~UCNField();
  
  virtual void GetFieldValue( const  double Point[4],
			      double *fieldArr ) const;
  virtual G4bool DoesFieldChangeEnergy() const{return true;}

private:
  
  std::vector<TField*> fields;
};

#endif
