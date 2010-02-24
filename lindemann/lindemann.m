function crosscorr = lindemann(insig,w_f,c_s,fs)
% LINDEMANN Calculates a binaural activation pattern
%   Usage: outsig = lindemann(insig,w_f,c_s,fs)
%
%   Input parameters:
%       insig       - binaural signal for which the cross-correlation
%                     should be calculated
%       w_f         - monaural sensitivity at the end of the delay line, 
%                     0 <= w_f < 1
%       c_s         - constant inhibition factor, 0 <= c_s <= 1 
%                     (0.0 = no inhibition)
%       fs          - sampling rate
%
%   Output parameters:
%       crosscorr   - A matrix containing the cross-correlation signal
%                     for every frequency channel fc and every time step n. 
%                     The format of this matrix is output(n,m,fc), where m
%                     denotes the correlation (delay line) time step.
%
%   LINDEMANN(insig,w_f,c_s,fs) calculates a binaural activity map for the
%   given insig using a cross-correlation (delay-line) mechanism. The
%   calculation is done for every frequency band in the range 5-40 Erb.
%
%   The steps of the binaural model to calculate the result are the
%   following:
%
%     1) The given stimulus is filtered using an erb bank to
%        get 36 frequency bands containing a stimulus waveform.
%
%     2) In a second step the auditory nerve is siumulated by extracting the
%        envelpoe using a first order low pass filter with a cutoff frequency
%        of 800 Hz and half-wave rectification.
%
%     3) Calculation of the cross-correlation between the left and right
%        channel.  This is done using the model described in Lindemann
%        (1986a). These are extensions to the delay line model of Jeffres (1948).
%        Lindemann has extended the delay line model of Jeffres (1948) by a
%        contralateral inhibition, which introduce the ILD to the model.  Also
%        monaural detectors were extended, to handle monaural signals (and some
%        stimuli with a split off of the lateralization image).  A detailed
%        description of these cross-correlation steps is given in the
%        binaural_correlation function.
%
%   See also: binncorr, plotlindemann, gammatone, filterbank
%
%   Demos: demo_lindemann
%
%R  lindemann1986a lindemann1986b gaik1993combined jeffres1948 hess2007phd

%   A first implementation of the Lindemann model in Matlab was done from
%   Wolfgang Hess and inspired this work.

% AUTHOR: Hagen Wierstorf


%% ------ Checking of input  parameters ---------------------------------

error(nargchk(4,4,nargin));

if ~isnumeric(insig)
    error('%s: insig has to be numeric!',upper(mfilename));
end


%% ------ Variables -----------------------------------------------------
% Highest and lowest frequency to use for the erbfilterbank (this gives us 
% 36 frequency channels, channel 5-40)
flow = erbtofreq(5);
fhigh = erbtofreq(40); 


% ------ Erb Bank -------------------------------------------------------
% Generate an erb filterbank for simulation of the frequncy -> place
% transformation of the cochlea. This generates erb filterbank coefficients
% with a range from flow to fhigh.
% NOTE: Lindemann uses a bandpass filterbank after Duifhuis (1972) and
% Blauert and Cobben (1978).
[b,a] = gammatone(erbspacebw(flow,fhigh),fs);
% Applying the erb filterbank to the signal
inoutsig = filterbank(b,a,insig);


%% ------ Cross-correlation computation ---------------------------------

% Extract the envelope, apply a half-wave rectification and calculate a
% running cross-correlation for every given frequency band
	
% ------ Haircell simulation -------
% Half-wave rectification and envelope extraction
inoutsig = ihcenvelope(inoutsig,fs,'lindemann');

% ------ cross-correlation ------
% Calculate the cross-correlation after Lindemann (1986a).
crosscorr = bincorr(inoutsig,w_f,c_s,fs);


