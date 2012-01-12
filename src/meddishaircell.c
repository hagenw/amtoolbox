#include "math.h"      

/* meddishaircell.c --  Meddis haircell model for monaural and
        binaural models (peripheral preprocessing).
        This file can be compiled by the matlab
        mex compiler. (c) 2001 Jeroen Breebaart
*/


      
void meddishaircell(
   double theVector[],    /* vector with data from VECTOR_IN */
   int fs,                /* samplerate of signal */
   int length,            /* length of vector */
   double theResult[]     /* vector with data for VECTOR_OUT */
   )
      
{

   int ii;

/* Parameters from Meddis' April 1990 JASA paper: */
   const double M=1;
   const double A=5;
   const double B=300;
   const double g=2000;
   const double y=5.05;
   const double l=2500;
   const double r=6580;
   const double x=66.31;
   const double h=50000;

/* internal variables */

   double kt;
   double spont;
   double q;
   double w;
   double temp;
   double replenish;
   double eject;
   double reuptake;
   double reprocess;
   double loss;

   double srate;

   printf("%i\n",fs);

   srate = 1.0*fs;

   kt = g*A/(A+B);
   spont = M*y*kt/(l*kt+y*(l+r));
   q=spont*(l+r)/kt;
   w=spont*r/x;

    
   /* do the MEDDISHAIRCELL thing! */
   for (ii=0; ii < length; ii++)
   {
      temp=(theVector[ii]+A+fabs(A+theVector[ii]))/2;
      kt=(g/srate)*temp/(temp+B);
      replenish=((y/srate)*(M-q)+fabs((y/srate)*(M-q)))/2;
      eject=kt*q;
      loss=(l/srate)*spont;
      reuptake=(r/srate)*spont;
      reprocess=(x/srate)*w;

      q=q+replenish-eject+reprocess;
      spont=spont+eject-loss-reuptake;
      w=w+reuptake-reprocess;
      theResult[ii]=spont*h;
    
   }
}
      
