function [crosscorr,t] = lindemann(insig,fs,varargin)
% LINDEMANN Calculates a binaural activation pattern
%   Usage: [crosscorr,t] = lindemann(insig,fs,c_s,w_f,M_f,T_int,N_1)
%          [crosscorr,t] = lindemann(insig,fs,c_s,w_f,M_f,T_int)
%          [crosscorr,t] = lindemann(insig,fs,c_s,w_f,M_f)
%          [crosscorr,t] = lindemann(insig,fs,c_s,w_f)
%          [crosscorr,t] = lindemann(insig,fs,c_s)
%          [crosscorr,t] = lindemann(insig,fs)
%
%   Input parameters:
%       insig       - binaural signal for which the cross-correlation
%                     should be calculated
%       fs          - sampling rate (Hz)
%       'stationary'- will set the default values of N_1=17640 and T_int=Inf, use
%                     this for stationary input signals
%       'dynamic'   - will set the default values of N_1=1 and T_int=5, use this
%                     for non-stationary input signals. This is the default.
%       c_s         - stationary inhibition factor, 0 <= c_s <= 1 
%                     (0.0 = no inhibition). Default: 0.3
%       w_f         - monaural sensitivity at the end of the delay line, 
%                     0 <= w_f < 1. Default: 0.035
%       M_f         - determines the decrease of the monaural sensitivity along 
%                     the delay line. Default: 6
%       T_int       - integration time window (ms). This is the memory of the 
%                     correlation process with exp(-1/T_int). Also this
%                     determines the time steps in the binaural activity map,
%                     because every time step T_int a new running
%                     cross-correlation is started, so every T_int we have a new
%                     result in crosscorr. You can set T_int = inf if you like
%                     to have no memory effects, then you will get only one
%                     time step in crosscorr. Default: T_int = 5~ms, T_int = Inf
%                     for 'stationary' (see above)
%       N_1         - Sample at which the first running cross-correlation should
%                     be started to avoid onset effects (see lindemann1986a p.
%                     1614). Default: 1, 17640 (200/f*fs with f=500 and fs=44100) 
%                     for 'stationary' (see above)
%
%   Output parameters:
%       crosscorr   - A matrix containing the cross-correlation signal
%                     for every frequency channel fc and every time step n. 
%                     The format of this matrix is output(n,m,fc), where m
%                     denotes the correlation (delay line) time step.
%       t           - time axis for the time steps n in crosscorr
%
%   LINDEMANN(insig,fs,c_s,w_f,M_f,T_int,N_1) calculates a binaural activity map
%   for the given insig using a cross-correlation (delay-line) mechanism. The
%   calculation is done for every frequency band in the range 5-40 Erb.
%
%   The steps of the binaural model to calculate the result are the
%   following:
%
%     1) The given stimulus is filtered using an erb bank to
%        get 36 frequency bands containing a stimulus waveform.
%
%     2) In a second step the auditory nerve is siumulated by extracting the
%        envelope using a first order low pass filter with a cutoff frequency
%        of 800 Hz and half-wave rectification.
%
%     3) Calculation of the cross-correlation between the left and right
%        channel.  This is done using the model described in Lindemann
%        (1986a) and Hess (2007). These are extensions to the delay line model 
%        of Jeffres (1948).
%
%        Lindemann has extended the delay line model of Jeffres (1948) by a
%        contralateral inhibition, which introduce the ILD to the model.  Also
%        monaural detectors were extended, to handle monaural signals (and some
%        stimuli with a split off of the lateralization image). Hess has
%        extented the output from the lindemann model to a binaural activity map
%        dependend on time, by using a running cross-correlation function.
%        This has been done here by starting a new running cross-correlation 
%        every time step T_int.  A detailed description of these cross-
%        correlation steps is given in the bincorr function.
%
%   See also: bincorr, plotlindemann, gammatone, filterbankz
%
%   Demos: demo_lindemann
%
%R  lindemann1986a lindemann1986b gaik1993combined jeffres1948 hess2007phd

%   A first implementation of the Lindemann model in Matlab was done from
%   Wolfgang Hess and inspired this work.

% AUTHOR: Hagen Wierstorf


%% ------ Checking of input  parameters ---------------------------------

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isnumeric(insig) || min(size(insig))~=2
    error('%s: insig has to be a numeric two channel signal!',upper(mfilename));
end

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
    error('%s: fs has to be a positive scalar!',upper(mfilename));
end

% For default values see lindemann1986a page 1613
% NOTE: I modified the default value for T_int from 10 to 5.

defnopos.flags.modus={'dynamic','stationary'};

defnopos.keyvals.c_s   = 0.3;
defnopos.keyvals.w_f   = 0.035;
defnopos.keyvals.M_f   = 6;
defnopos.keyvals.T_int = 5;
defnopos.keyvals.N_1   = 1;

% Parse the command line
[flags,keyvals,c_s,w_f,M_f,T_int,N_1]  = ...
    ltfatarghelper({'c_s','w_f','M_f','T_int','N_1'},defnopos,varargin);

if flags.do_stationary
  T_int = Inf;
  N_1   = 17640;  
end;

if ( ~isnumeric(c_s) || ~isscalar(c_s) || c_s<0 || c_s>1 )
    error('%s: 0 <= c_s <= 1, but c_s = %.1f',upper(mfilename),c_s);
end

if ( ~isnumeric(w_f) || ~isscalar(w_f) || w_f<0 || w_f>=1 )
    error('%s: 0 <= w_f < 1, but w_f = %.1f',upper(mfilename),w_f);
end

if ( ~isnumeric(M_f) || ~isscalar(M_f) || M_f<=0 )
    error('%s: M_f has to be a positive scalar!',upper(mfilename));
end

if ( ~isnumeric(T_int) || ~isscalar(T_int) || T_int<=0 )
    error('%s: T_int has to be a positive scalar!',upper(mfilename));
end

if ( ~isnumeric(N_1) || ~isscalar(N_1) || N_1<=0 )
    error('%s: N_1 has to be a positive scalar!',upper(mfilename));
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
%
% FIXME: There is currently an error in the gammatone function for real valued
% filters, so use complex valued filters instead.
[b,a] = gammatone(erbspacebw(flow,fhigh),fs,'complex');
% Applying the erb filterbank to the signal
inoutsig = 2*real(filterbankz(b,a,insig));


%% ------ Cross-correlation computation ---------------------------------

% Extract the envelope, apply a half-wave rectification and calculate a
% running cross-correlation for every given frequency band
	
% ------ Haircell simulation -------
% Half-wave rectification and envelope extraction
inoutsig = ihcenvelope(inoutsig,fs,'lindemann');

% ------ Cross-correlation ------
% Calculate the cross-correlation after Lindemann (1986a).
[crosscorr,t] = bincorr(inoutsig,fs,c_s,w_f,M_f,T_int,N_1);


