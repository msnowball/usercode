#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "HiggsCSandWidthSM4.cc"

using namespace std;

int main()
{

  ofstream fileOut;
  char* fileName_[6] = {"Total_cs_pdf.txt","GluGlu_cs_pdf.txt","VBF_cs_pdf.txt","WH_cs_pdf.txt","ZH_cs_pdf.txt","ttH_cs_pdf.txt"};

  HiggsCSandWidthSM4 *myCSW = new HiggsCSandWidthSM4();

  for( int i = 0; i < 6; i++)
    {
      
      fileOut.open(fileName_[i]);
      fileOut << " mH      CS   CSErr Plus%   CSErr Plus   CSErr Minus%  CSErr Minus " << endl;

      double mH;
      double CS;
      double CSErrPlusPercent, CSErrMinusPercent;
      double CSErrPlus, CSErrMinus;
      double sqrts = 7;
      bool spline = true;

      for( double j = 110; j < 140; j += 0.5)
	{

	  mH = j;
	  CS = myCSW->HiggsCS(i,mH,sqrts,spline);
	  CSErrPlusPercent = myCSW->HiggsCSpdfErrPlus(i,mH,sqrts)*100;
	  CSErrPlus = CS+myCSW->HiggsCSpdfErrPlus(i,mH,sqrts)*CS;
	  CSErrMinusPercent = myCSW->HiggsCSpdfErrMinus(i,mH,sqrts)*100;
	  CSErrMinus = CS+myCSW->HiggsCSpdfErrMinus(i,mH,sqrts)*CS;
	  
	  fileOut.width(6);
	  fileOut << mH << "  " << CS << "  " << CSErrPlusPercent << "      " << CSErrPlus << "        " << CSErrMinusPercent << "       " << CSErrMinus << endl;

	}



      for( double k = 140; k < 160; k++ )
	{

          mH = k;
          CS = myCSW->HiggsCS(i,mH,sqrts,spline);
          CSErrPlusPercent = myCSW->HiggsCSpdfErrPlus(i,mH,sqrts)*100;
          CSErrPlus = CS+myCSW->HiggsCSpdfErrPlus(i,mH,sqrts)*CS;
          CSErrMinusPercent = myCSW->HiggsCSpdfErrMinus(i,mH,sqrts)*100;
          CSErrMinus = CS+myCSW->HiggsCSpdfErrMinus(i,mH,sqrts)*CS;
	  
	  fileOut.width(6);
	  fileOut << mH << "  " << CS << "  " << CSErrPlusPercent << "      " << CSErrPlus << "        " << CSErrMinusPercent << "       " << CSErrMinus << endl;


	}

      for( double l = 160; l < 290; l += 2)
	{
          mH = l;
          CS = myCSW->HiggsCS(i,mH,sqrts,spline);
          CSErrPlusPercent = myCSW->HiggsCSpdfErrPlus(i,mH,sqrts)*100;
          CSErrPlus = CS+myCSW->HiggsCSpdfErrPlus(i,mH,sqrts)*CS;
          CSErrMinusPercent = myCSW->HiggsCSpdfErrMinus(i,mH,sqrts)*100;
          CSErrMinus = CS+myCSW->HiggsCSpdfErrMinus(i,mH,sqrts)*CS;

	  fileOut.width(6);
          fileOut << mH << "  " << CS << "  " << CSErrPlusPercent << "      " << CSErrPlus << "        " << CSErrMinusPercent << "       " << CSErrMinus << endl;


	}

      for( double m = 290; m < 350; m += 5)
	{

          mH = m;
          CS = myCSW->HiggsCS(i,mH,sqrts,spline);
          CSErrPlusPercent = myCSW->HiggsCSpdfErrPlus(i,mH,sqrts)*100;
          CSErrPlus = CS+myCSW->HiggsCSpdfErrPlus(i,mH,sqrts)*CS;
          CSErrMinusPercent = myCSW->HiggsCSpdfErrMinus(i,mH,sqrts)*100;
          CSErrMinus = CS+myCSW->HiggsCSpdfErrMinus(i,mH,sqrts)*CS;

	  fileOut.width(6);
          fileOut << mH << "  " << CS << "  " << CSErrPlusPercent << "      " << CSErrPlus << "        " << CSErrMinusPercent << "       " << CSErrMinus << endl;


	}

      for( double n = 350; n < 400; n += 10)
	{

          mH = n;
          CS = myCSW->HiggsCS(i,mH,sqrts,spline);
          CSErrPlusPercent = myCSW->HiggsCSpdfErrPlus(i,mH,sqrts)*100;
          CSErrPlus = CS+myCSW->HiggsCSpdfErrPlus(i,mH,sqrts)*CS;
          CSErrMinusPercent = myCSW->HiggsCSpdfErrMinus(i,mH,sqrts)*100;
          CSErrMinus = CS+myCSW->HiggsCSpdfErrMinus(i,mH,sqrts)*CS;

	  fileOut.width(6);
          fileOut << mH << "  " << CS << "  " << CSErrPlusPercent << "      " << CSErrPlus << "        " << CSErrMinusPercent << "       " << CSErrMinus << endl;


	}

      for( double q = 400; q <= 1000; q += 20 )
	{

          mH = q;
          CS = myCSW->HiggsCS(i,mH,sqrts,spline);
          CSErrPlusPercent = myCSW->HiggsCSpdfErrPlus(i,mH,sqrts)*100;
          CSErrPlus = CS+myCSW->HiggsCSpdfErrPlus(i,mH,sqrts)*CS;
          CSErrMinusPercent = myCSW->HiggsCSpdfErrMinus(i,mH,sqrts)*100;
          CSErrMinus = CS+myCSW->HiggsCSpdfErrMinus(i,mH,sqrts)*CS;

	  fileOut.width(6);
          fileOut << mH << "  " << CS << "  " << CSErrPlusPercent << "      " << CSErrPlus << "        " << CSErrMinusPercent << "       " << CSErrMinus << endl;


	}

      fileOut.close();

    }




  return 0;

}
