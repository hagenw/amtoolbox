function [crosscorr,t,ild,cfreq] = lindemann(insig,fs,varargin)
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
%
%   Output parameters:
%       crosscorr   - A matrix containing the cross-correlation signal
%                     for every frequency channel fc and every time step n.
%                     The format of this matrix is output(n,m,fc), where m
%                     denotes the correlation (delay line) time step.
%       t           - time axis for the time steps n in crosscorr
%       ild         - interaural level difference (ILD) for every freqeuncy
%                     channel fc
%       cfreq       - center frequencies of every frequency channel
%
%   LINDEMANN(insig,fs) calculates a binaural activity map for the given
%   insig using a cross-correlation (delay-line) mechanism. The calculation
%   is done for every frequency band in the range 5-40 Erb.
%
%   Lindemann has extended the delay line model of Jeffres (1948) by a
%   contralateral inhibition, which introduce the ILD to the model.  Also
%   monaural detectors were extended, to handle monaural signals (and some
%   stimuli with a split off of the lateralization image). Hess has
%   extented the output from the lindemann model to a binaural activity map
%   dependend on time, by using a running cross-correlation function.
%   This has been done here by starting a new running cross-correlation
%   every time step T_int.  A detailed description of these cross-
%   correlation steps is given in the lindemannbincorr function.

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
%   You may supply any flags or key/value pairs of the AUDIORYFILTERBANK,
%   IHCENVELOPE or LINDEMANNBINCORR at the end of the line of input
%   arguments.
%
%   See also: lindemannbincorr, plotlindemann, gammatone, ufilterbankz
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

% Parse the command line and load default parameters
% For default values see lindemann1986a page 1613
% NOTE: I modified the default value for T_int from 10 to 5.
definput.import={'auditoryfilterbank','ihcenvelope','lindemannbincorr'};
% Highest and lowest frequency to use for the erbfilterbank (this gives us
% 36 frequency channels, channel 5-40)
definput.importdefaults = ...
    {'flow',erbtofreq(5),'fhigh',erbtofreq(40),'ihc_lindemann'};
[flags,keyvals,c_s,w_f,M_f,T_int,N_1]  = ...
    ltfatarghelper({'c_s','w_f','M_f','T_int','N_1'},definput,varargin);


%% ------ Computation ---------------------------------------------------

% ------ Erb Bank -------------------------------------------------------
% Apply the auditory filterbank
% NOTE: Lindemann uses a bandpass filterbank after Duifhuis (1972) and
% Blauert and Cobben (1978).
[inoutsig,cfreq] = auditoryfilterbank(insig,fs,'argimport',flags,keyvals);

% ------ ILD ------------------------------------------------------------
% Calculate the interaural level difference (ILD) for every frequency channel
% NOTE: this was not part of the original Lindemann model
ild = interauralleveldifference(inoutsig(:,:,1),inoutsig(:,:,2));

% ------ Cross-correlation computation ---------------------------------
% Extract the envelope, apply a half-wave rectification and calculate a
% running cross-correlation for every given frequency band
% ------ Haircell simulation -------
% Half-wave rectification and envelope extraction
inoutsig = ihcenvelope(inoutsig,fs,'argimport',flags,keyvals);
% ------ Cross-correlation ------
% Calculate the cross-correlation after Lindemann (1986a).
[crosscorr,t] = lindemannbincorr(inoutsig,fs,c_s,w_f,M_f,T_int,N_1);
