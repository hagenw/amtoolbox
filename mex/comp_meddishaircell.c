#include "mex.h"

#include "../src/meddishaircell.c"


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
      
 
void mexFunction(
   int     nlhs,
   mxArray  *plhs[],
   int     nrhs,
   const mxArray  *prhs[]
   )
{

   int  fs;
   int  m,n;            /* temporary array size holders */


   m=mxGetM(prhs[0]);  /* Assigning result of mxGetM and
			    mxGetN to variables */
   n=mxGetN(prhs[0]);  /* with simpler names makes for
			    easier code reading.  */
   
/*   Step 1: Error Checking Step 1a: is nlhs 1?  If not,
     generate an error message and exit mexample (mexErrMsgTxt
     does this for us!) */
   if (nlhs!=1)
      mexErrMsgTxt("meddishaircell requires one output argument.");
      
   /*   Step 1b: is nrhs 2? */
   if (nrhs!=2)
      mexErrMsgTxt("meddishaircell requires two input arguments, \
    prhs[0] and SAMPLERATE");
   
   /* step 1c: is prhs[0] a x by 1 row? */
   if (m<2)
      mexErrMsgTxt("meddishaircell requires one column as input signal.");
    
   if (n!=1)
      mexErrMsgTxt("meddishaircell requires one column as input signal.");


   fs=(int)mxGetScalar(prhs[1]);
  
   plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);
      
   meddishaircell(mxGetPr(prhs[0]), fs, m, 
		  mxGetPr(plhs[0]));
 
}
