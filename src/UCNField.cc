#include "G4SystemOfUnits.hh"

#include "G4Field.hh"
#include "UCNGlobals.hh"
#include "UCN2DField.hh"
#include "UCN3DField.hh"
#include "UCNConductorField.hh"
#include "UCNField.hh"


UCNField::UCNField(std::vector<TField*> FIELDS)
{
  fields=FIELDS;
}

UCNField::~UCNField()
{;}

void UCNField::GetFieldValue(const double Point[4],double *fieldArr) const
{
  double x, y, z, t;
  double V;
  double Ei[3];
  double B[4][4];

  x=Point[0]/m;
  y=Point[1]/m;
  z=Point[2]/m;
  t=Point[3]/s;

  for (int k = 0; k < 4; k++)
    for (int j = 0; j < 4; j++)
      B[k][j] = 0;
  Ei[0] = Ei[1] = Ei[2] = V = 0;

  //for (std::vector<TField*>::iterator i = fields.begin(); i != fields.end(); i++){
  for (std::vector<TField*>::const_iterator i = fields.begin(); i != fields.end(); i++){
    (*i)->BField(x, y, z, t, B);
    (*i)->EField(x, y, z, t, V, Ei);
  }
  fieldArr[0]=B[0][0]*tesla;
  fieldArr[1]=B[1][0]*tesla;
  fieldArr[2]=B[2][0]*tesla;
  fieldArr[3]=Ei[0]*volt/m;
  fieldArr[4]=Ei[1]*volt/m;
  fieldArr[5]=Ei[2]*volt/m;
}

