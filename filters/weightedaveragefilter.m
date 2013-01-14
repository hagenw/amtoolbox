function outsig = weightedaveragefilter(insig,weight,fs,timeconst)
%WEIGHTEDAVERAGEFILTER Compute the weighted or self-weighted average
%                 
%   Usage: outsig = weightedaveragefilter(insig,weight,fs,timeconst)
%
%   Input parameters:
%        insig     : signal from which the average is computed. Optionally,
%                    can be the same as the weight resulting in
%                    self-weighted average
%        weight    : signal to be used as weight in the computation
%        fs        : sampling rate
%        timeconst : time constant specifying the first-order IIR filter
%
%   Output parameters:
%        outsig    : output signal
%
%   This function computes the either the weighted or the self-weighted
%   average of the input signal using a first-order IIR filter whose time
%   constant is specified by the *timeconst* argument. More details about
%   the conputation can be found in Pulkki, Hirvonen 2009 (Sec. 3.2.3)
%
%   See also: takanen2013mso, takanen2013lso, takanen2013wbmso,
%             takanen2013onsetenhancement, takanen2013contracomparison
%
%   References: pulkki2009 takanen2013a

%   AUTHOR: Marko Takanen, Olli Santala, Ville Pulkki
%
%   COPYRIGHT (C) 2013 Aalto University
%                      School of Electrical Engineering
%                      Department of Signal Processing and Acoustics
%                      Espoo, Finland


%% ------ Computation ----------------------------------------------------

%setting the parameters for a first-order IIR filter
B = 1-exp(-1/(fs*timeconst));
A = [1 -exp(-1/(fs*timeconst))];

%filter the weighted and self-weighted signals 
weighted = filter(B,A,(insig.*(weight.^2)));
selfWeighted = filter(B,A,(weight.^2));

%derive the output
outsig = weighted./(selfWeighted+1e-30);