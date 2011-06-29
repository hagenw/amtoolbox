/* fadapt.c --  IHC model for monaural and
        binaural models (peripheral preprocessing).
        This file can be compiled by the matlab
        mex compiler. (c) 2000 Jeroen Breebaart
*/


      
#include "mex.h"
#include "math.h"      
      
void workFcn(
   double srate[],    /* samplerate of signal */
   int length,        /* length of vector */
   double theVector[],    /* vector with data from VECTOR_IN */
   double theResult[]     /* vector with data for VECTOR_OUT */
        )
      
{
    int teller,order;
    const double pi = 3.1415926535898;    
    double b;
    
    
    /* do the IHC thing! */
    for (teller=0; teller < length; teller++)
        {theResult[teller] = theVector[teller];
        if (theResult[teller] < 0){theResult[teller]=0;} 
        }

    /* lowpass filtering:*/
    b = 2 - cos(2*pi*2000 / srate[0]) - pow((pow((cos(2*pi*2000 / srate[0])-2),2)-1),0.5);
    order = 5;
    while (order > 0)
        {
        theResult[0] = theResult[0]*(1-b);
        for (teller=1; teller < length; teller++)
            {theResult[teller] = theResult[teller] * (1-b) + b * theResult[teller-1];
            }
        order = order - 1;
        }
    




}
      
/*   Now, define the gateway function, i.e., mexFunction.Below
     is the standard, predeclared header to mexFunction.  nlhs
     and nrhs are the number of left-hand and right-hand side
     arguments that mexample was called with from within MATLAB.
     In this example, nlhs equals 1 and nrhs should equal 2.  If
     not, then the user has called mexample the wrong way and
     should be informed of this.  plhs and prhs are arrays which
     contain the pointers to the MATLAB arrays, which are
     stored in a C struct called an Array.  prhs is an array of
     length rhs,and its pointers point to valid input data.
     plhs is an array of length nlhs, and its pointers point to
     invalid data (i.e., garbage).  It is the job of mexFunction
     to fill plhs with valid data.
      
     First, define the following values.  This makes it much
     easier to change the order of inputs to mexample, should we
     want to change the function later.  In addition, it makes
     the code easier to read. */
      

#define VECTOR_IN prhs[0]
#define VECTOR_OUT plhs[0]
#define fs prhs[1]
 
void mexFunction(
        int     nlhs,
        mxArray  *plhs[],
        int     nrhs,
        const mxArray  *prhs[]
        )
{


int    m,n;            /* temporary array size holders */



m=mxGetM(VECTOR_IN);  /* Assigning result of mxGetM and
                           mxGetN to variables */
n=mxGetN(VECTOR_IN);  /* with simpler names makes for
                           easier code reading.  */
   
/*   Step 1: Error Checking Step 1a: is nlhs 1?  If not,
     generate an error message and exit mexample (mexErrMsgTxt
     does this for us!) */
  if (nlhs!=1)
    mexErrMsgTxt("ihc requires one output argument.");
      
  /*   Step 1b: is nrhs 2? */
  if (nrhs!=2)
    mexErrMsgTxt("ihc requires two input arguments, \
    VECTOR_IN and SAMPLERATE");
   
   /* step 1c: is VECTOR_IN a x by 1 row? */
  if (m<2)
    mexErrMsgTxt("ihc requires one column as input signal.");
    
  if (n!=1)
    mexErrMsgTxt("ihc requires one column as input signal.");

   
  /*   Step 2:  Allocate memory for return argument(s) */
  VECTOR_OUT = mxCreateDoubleMatrix(m, 1, mxREAL);  /* create a mx1
                                                       full, numeric,
                                                       real-valued
                                                       array */
      

  /*   mxGetPr returns a pointer to the real part of the array
       ARRAY_IN.  In the line above, it is treated as the one-
       dimensional array mentioned in the previous comment.  */
      
  /*   Step 4:  Call workFcn function */
  workFcn(mxGetPr(fs), m, mxGetPr(VECTOR_IN), mxGetPr(VECTOR_OUT));
 
}
