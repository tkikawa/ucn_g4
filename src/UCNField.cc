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

  x = Point[0]/m;
  y = Point[1]/m;
  z = Point[2]/m;
  t = Point[3]/s;

  for (int k = 0; k < 4; k++)
    for (int j = 0; j < 4; j++)
      B[k][j] = 0;
  Ei[0] = Ei[1] = Ei[2] = V = 0;
  
  for (std::vector<TField*>::const_iterator i = fields.begin(); i != fields.end(); i++){
    (*i)->BField(x, y, z, t, B);
    (*i)->EField(x, y, z, t, V, Ei);
  }

  fieldArr[0] = B[0][0]*tesla;  //Bx
  fieldArr[1] = B[1][0]*tesla;  //By
  fieldArr[2] = B[2][0]*tesla;  //Bz
  fieldArr[3] = Ei[0]*volt/m;   //Ex
  fieldArr[4] = Ei[1]*volt/m;   //Ey
  fieldArr[5] = Ei[2]*volt/m;   //Ez
  fieldArr[6] = 0;              //gx
  fieldArr[7] = 0;              //gy
  fieldArr[8] = -9.81*m/s/s;    //gz
  //fieldArr[8] = 0;    //gz
  fieldArr[9] = B[0][1]*tesla/m;//dBx/dx
  fieldArr[10]= B[1][1]*tesla/m;//dBy/dx
  fieldArr[11]= B[2][1]*tesla/m;//dBz/dx
  fieldArr[12]= B[0][2]*tesla/m;//dBx/dy
  fieldArr[13]= B[1][2]*tesla/m;//dBy/dy
  fieldArr[14]= B[2][2]*tesla/m;//dBz/dy
  fieldArr[15]= B[0][3]*tesla/m;//dBx/dz
  fieldArr[16]= B[1][3]*tesla/m;//dBy/dz
  fieldArr[17]= B[2][3]*tesla/m;//dBz/dz
}

void UCNField::GetCurrentFieldValue(const double t, const double x, const double y, const double z, double B[4][4], double *Ei, double &V)
{
  for (int k = 0; k < 4; k++)
    for (int j = 0; j < 4; j++)
      B[k][j] = 0;
  Ei[0] = Ei[1] = Ei[2] = V = 0;

  for (std::vector<TField*>::const_iterator i = fields.begin(); i != fields.end(); i++){
    (*i)->BField(x, y, z, t, B);
    (*i)->EField(x, y, z, t, V, Ei);
  }



  B[3][0] = sqrt(B[0][0]*B[0][0] + B[1][0]*B[1][0] + B[2][0]*B[2][0]); // absolute value of B-Vector
  if (B[3][0]>1e-31)
    {
      B[3][1] = (B[0][0]*B[0][1] + B[1][0]*B[1][1] + B[2][0]*B[2][1])/B[3][0]; // derivatives of absolute value
      B[3][2] = (B[0][0]*B[0][2] + B[1][0]*B[1][2] + B[2][0]*B[2][2])/B[3][0];
      B[3][3] = (B[0][0]*B[0][3] + B[1][0]*B[1][3] + B[2][0]*B[2][3])/B[3][0];
    }

  

}
