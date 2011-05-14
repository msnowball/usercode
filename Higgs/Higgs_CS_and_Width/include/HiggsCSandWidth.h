#ifndef HIGGSCSANDWIDTH_H
#define HIGGSCSANDWIDTH_H

#define PI 3.14159

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>


/**********************************************************/
/*            Class for Higgs Width and CS                */
/*                                                        */
/*  All numbers for CS and width are taken from official  */
/*  numbers on Higgs CS Twiki (Spring 2011)               */
/*                                                        */
/*  These numbers are taken into memory and a simple      */
/*  linear interpolation is done.                         */
/*                                                        */
/*  For any invalid process or mH out of range, -1 will   */
/*  be returned.                                          */
/*                                                        */
/*    Written by:                                         */
/*         Matt Snowball                                  */
/*         University of Florida                          */
/*         snowball@phys.ufl.edu                          */
/**********************************************************/



class HiggsCSandWidth
{

 public:

  HiggsCSandWidth();
  ~HiggsCSandWidth();

  double HiggsCS(int ID, double mH, double sqrts);
  double HiggsCSErrPlus(int ID, double mH, double sqrts);
  double HiggsCSErrMinus(int ID, double mH, double sqrts);
  double HiggsCSscaleErrPlus(int ID, double mH, double sqrts);
  double HiggsCSscaleErrMinus(int ID, double mH, double sqrts);
  double HiggsCSpdfErrPlus(int ID, double mH, double sqrts);
  double HiggsCSpdfErrMinus(int ID, double mH, double sqrts);
  double HiggsWidth(int ID,double mH);


 private:

  double scratchMass;
  double BR[26][103];
  double CS[6][50];
  double CSerrPlus[6][50];
  double CSerrMinus[6][50];
  double CSscaleErrPlus[6][50];
  double CSscaleErrMinus[6][50];
  double CSpdfErrPlus[6][50];
  double CSpdfErrMinus[6][50];

  std::string FileLoc;


};

#endif
