/* prgtf.c - C-code for a near-Perfect Reconstruction GammaTone Filterbank (PRGTF)
   Programmed by Jeroen Breebaart, (C) 2003. This code implements a 3th-order, <=40 chan
   gammatone filterbank with filters at ERB spacing (Glasberg & Moore estimates).
   The current filterbank design includes a delay of 1024 samples.

   To start the filterbank, first call the initializer: initGtf();
   A new set of output samples is put in the array Out_real and Out_imag whenever
   nextCycle(input sample) is called.

   Summing the real part of the outputs will result in near-perfect reconstruction of the
   original input signal (within 1 dB over the complete frequency range).
   The absolute value of the (complex) output gives the Hilbert envelope of the filtered signals.
*/

   


# include "math.h"
# include "mex.h"


/* Define global variables: */

int fs;				                	/* Sample rate [Hz] */
double T;				            	/* Overall group delay [s] */
const n = 3;						/* order of gammatone filters [-]*/
double cf[40];						/* Center frequencies of channels [Hz] */
double b[40];						/* bandwidth parameter of channels [rad/s] */
double Delta[40];					/* Delay of maximum envelope [s] */
double Theta[40];					/* necessary phase delay [rad] */
double State1_re[40];					/* states of first LP filter, complex */
double State1_im[40];
double State2_re[40];					/* states of second LP filter, complex */
double State2_im[40];
double State3_re[40];					/* states of third LP filter, complex */
double State3_im[40];
double Out_real[40];			    		/* successive output samples for each filter, real part */
double Out_imag[40];                			/* successive output samples for each filter, imaginary part */
double Xponent[40];					/* exp(b/fs) */
int Delay[40];						/* delay to apply to filters [samples] */
double snum;						/* sample number of input */
double sampbuf[1024];					/* sample buffer for delays */
int bpoint;						/* pointer where to write data to buffer */
const double pi = 3.1415926535898;			/* if you don't understand these, */
const double e = 2.718281828459;			/* you should not read this code... */
const double sqrt2 = 1.41421356237310;			/* sqrt 2*/
int numfilt;                        			/* total number of filters used */
double spacing;                     			/* filter spacing in ERB rates */
double cosTable[1024];					/* table for fast cosine/sine computation */
const double cosConv = 162.9746617261008;		/* conversion from angle to cosine index */

/* mycos function */

double mycos(double x)
{
	int idx = (int)(x * cosConv);
	idx = (idx & 1023);
	return cosTable[idx];	
}

double mysin(double x)
{
	return mycos(pi/2 - x);
}

/* helper function for factors */

int factorial(int input)
{
	int output,teller;
	output=1;
	if (input > 1)
    {
		for (teller = 1; teller <= input; teller++){output = output * teller;}
    }

	return output;

}

/* stream one sample through filterbank */

void nextCycle(double input)
{
	int teller;
	double temp_re, temp_im;
	sampbuf[bpoint]=input;			/* point new input sample in buffer */
	bpoint=(bpoint+1)&1023;
	snum=snum+1;

	for (teller=0; teller < numfilt ; teller++)
	{
		/* read delayed sample from buffer and shift in frequency: */
		temp_re = sampbuf[(bpoint-Delay[teller])&1023] * mycos(-2*pi*cf[teller]*snum / fs + Theta[teller]);
        	temp_im = sampbuf[(bpoint-Delay[teller])&1023] * mysin(-2*pi*cf[teller]*snum / fs + Theta[teller]);
		
		/* lowpass filter: 3 cascaded 1st order IIR sections */
		temp_re = temp_re * (1-Xponent[teller]) + Xponent[teller] * State1_re[teller];
        	temp_im = temp_im * (1-Xponent[teller]) + Xponent[teller] * State1_im[teller];
		State1_re[teller]=temp_re;
		State1_im[teller]=temp_im;

		temp_re = temp_re * (1-Xponent[teller]) + Xponent[teller] * State2_re[teller];
        	temp_im = temp_im * (1-Xponent[teller]) + Xponent[teller] * State2_im[teller];
		State2_re[teller]=temp_re;
		State2_im[teller]=temp_im;

		temp_re = temp_re * (1-Xponent[teller]) + Xponent[teller] * State3_re[teller];
        	temp_im = temp_im * (1-Xponent[teller]) + Xponent[teller] * State3_im[teller];
		State3_re[teller]=temp_re;
		State3_im[teller]=temp_im;

		/* compute output */
        	Out_real[teller] = sqrt2 * temp_re * mycos(2*pi*cf[teller] * snum / fs) - sqrt(2) * temp_im * mysin(2*pi*cf[teller]*snum / fs);      
        	Out_imag[teller] = sqrt2 * temp_re * mysin(2*pi*cf[teller] * snum / fs) + sqrt(2) * temp_im * mycos(2*pi*cf[teller]*snum / fs);
	
	}

}

/* GTF initializer */

void initGtf()
{

	int teller;
	double erb;
    	T=1024/fs;
	snum = 0;

	/* initialize variables: */
	for (teller = 1; teller < numfilt+1; teller++)
	{
	    	erb=(teller-1)*spacing+1;
		cf[teller-1]=(pow(e,0.11*(erb))-1)/.00437;
		b[teller-1] = -2*spacing*pi*24.7*(pow(e,0.11*(erb)))*pow((factorial(n-1)),2)/(pi * factorial(2*n-2) * pow(2,-2*n+2));
		Delta[teller-1]=(-n+1)/b[teller-1];
		
		/* check if Delta is sufficiently small. Otherwise increase bandwidth of current filter to match maximum delay */
		if (Delta[teller-1]>1000)
		{
		    b[teller-1]=(-n+1)/1000*fs;
		    Delta[teller-1]=1000/fs;
		};
		
		Theta[teller-1]=-2*pi*Delta[teller-1]*cf[teller-1];
		Delay[teller-1]=floor((T-Delta[teller-1])*fs+0.5);
		
		
		State1_re[teller-1]=0;
		State2_re[teller-1]=0;
		State3_re[teller-1]=0;
		State1_im[teller-1]=0;
		State2_im[teller-1]=0;
		State3_im[teller-1]=0;
		Xponent[teller-1]=pow(e,b[teller-1]/fs);
		
	}
	
	

	for (teller=0; teller < 1024 ; teller++)
	{
		sampbuf[teller]=0;
	}

	bpoint=0;

	for (teller=0; teller < 1024 ; teller++)
	{
		cosTable[teller]=cos(2*pi*teller/1024);	
	}

}



/* work function: streams all input signals through filterbank */

void workFcn(
   const int length,             /* length of vector */
   double theVector[],           /* vector with real data from VECTOR_IN */
   double theResult_real[],	 /* real data from VECTOR_OUT */
   double theResult_imag[]       /* imag data from VECTOR_OUT */
)
{
	int sampnum, erbnum;

	initGtf();
    
	
	/* start streaming samples. The first 1024 samples are used to fill the buffer (delay compensation)*/
	
	for (sampnum = 0; sampnum < 1024 ; sampnum++)
	{
		nextCycle(theVector[sampnum]);
	}	
	
	
	for (sampnum = 1024; sampnum < length ; sampnum++)
	{
		nextCycle(theVector[sampnum]);
		for (erbnum=0 ; erbnum < numfilt ; erbnum++)
		{
			theResult_real[erbnum*length+sampnum-1024]=Out_real[erbnum];
			theResult_imag[erbnum*length+sampnum-1024]=Out_imag[erbnum];
		}
	}


	for (sampnum = length; sampnum < length+1024 ; sampnum++)
	{
		nextCycle(0);
		for (erbnum=0 ; erbnum < numfilt ; erbnum++)
		{
			theResult_real[erbnum*length+sampnum-1024]=Out_real[erbnum];
			theResult_imag[erbnum*length+sampnum-1024]=Out_imag[erbnum];
		}
	}
	
	
	
}





#define VECTOR_IN prhs[0]
#define VECTOR_OUT plhs[0]

 
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
    mexErrMsgTxt("prgtf requires one output argument.");
      
  /*   Step 1b: is nrhs 4? */
  if (nrhs!=4)
    mexErrMsgTxt("prgtf requires four input arguments:  VECTOR_IN, SAMPLERATE, NUMFILT, SPACING");
   
   /* step 1c: is VECTOR_IN a x by 1 row? */
  if (m<2)
    mexErrMsgTxt("prgtf requires one column as input signal.");
    
  if (n!=1)
    mexErrMsgTxt("prgtf requires one column as input signal.");


   /* step 2: read variables and to appropriate initialization */
   numfilt      = *mxGetPr(prhs[2]);
   fs           = *mxGetPr(prhs[1]);
   spacing      = *mxGetPr(prhs[3]);
   
   
   if( (pow(e,0.11*(numfilt*spacing))-1)/.00437    > fs/2)
   mexErrMsgTxt("Center frequency of highest filter is beyond Nyquist frequency. Lower spacing or number of filters.");
  
   if (numfilt > 40)
   mexErrMsgTxt("The maximum allowed number of filters amounts 40.");
   
  /*   Step 3:  Allocate memory for return argument(s) */
  VECTOR_OUT = mxCreateDoubleMatrix(m, numfilt, mxCOMPLEX);  /* create a mx1
                                                       full, numeric,
                                                       complex-valued
                                                       array */

       
       
  /*   Step 4:  Call workFcn function */
  workFcn(m, mxGetPr(VECTOR_IN), mxGetPr(VECTOR_OUT),mxGetPi(VECTOR_OUT));
 
}




