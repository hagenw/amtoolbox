function [ itdindices envindices ildindices ] = hrtf_analysis(insig,fs)
%__________________________________________________________________________
%
% hrtf_analysis( filtered_hrtf_signal )
%
% Calculates the ITD, the envelope correlation and the ILD for all bands 
% between the left and right channel of the given HRTF signal for a single 
% direction.
%
% Options:
%   filtered_hrtf_signal    - Struct containing the (ERB) band filtered
%                             HRTF signals:
%                                filtered_hrtf_signal.signal_left(t,band)
%                                filtered_hrtf_signal.signal_right(t,band)
%
% Output:
%   itd_indices             - Array containing the positions of the ITDs 
%                             (maxima of the cross-correlation function) 
%                             for every band
%   envelope_indices        - Array containing the positions of the maxima
%                             of the cross-correlation functions for the
%                             envelopes of the signals in every band
%   ild_indices             - Array containing the ILDs for every band 
%
%
% This function is mainly the implementation of a algorithm described in
% Gaik (1993). He calculates the ILD and ITD for measured HRTFs to identify
% naturally occuring combinations of ILD and ITD.
%
% Gaik (1993):
% "Combined evaluation of interaural time and intensity differences:
% Psychoacoustic results and computer modeling", Werner Gaik, 1993 JASA
%
% see also: filter_bank, spatialize
%__________________________________________________________________________
%bi.mo                                                                v.0.3

%
% Hagen Wierstorf
% T-Labs, Berlin
% hagen.wierstorf@telekom.de
% 2009/05/15
%

% Calculate the correlations for all (erb) bands
for band = 1:length(insig.signal_left(1,:))
    
    
		% ------ ITD correlation --------------------------------------------------
    
    sigl = insig.signal_left(:,band);
	  sigr = insig.signal_right(:,band);
    
    % Calculate cross-correlation between left and right channel
    %>(this is called itd_correlation, because the HRTFs are already
    %>applied and the left and right channel have an ITD)
    itdcorr(:,band) = xcorr(sigl,sigr,48);
    
    % Search for the maximum of the correlation in this band and store its 
    %>index
    [maximum,itdmaxindex(band)] = max(itdcorr(:,band));
    
    
		% ------ Envelope correlation ---------------------------------------------
    
    % Extract the envelopes of the right and the left signal using the
    %>hilbert transformation
    envl = calculate_envelope(sigl,0,fs,'m2');                                           
    envr = calculate_envelope(sigr,0,fs,'m2');
    
		% Calculate the cross-correlation between the left and right envelope
    envcorr(:,band) = xcorr(envl,envr,48);
    
    % Search for the maximum of the correlation in this band and store its
    %>index
    [maximum,envmaxindex(band)] = max(envcorr(:,band));
    
    
		% ------ ILD --------------------------------------------------------------
    
		% Calculate the interaural level difference (ILD) for this band.
		%>(use 10*log10, because we have the squares of the signals!)
    ild(band) = 10.*log10( mean(sigl.^2) ./ mean(sigr.^2) ); 
    
end % of for


%% Return the results
itdindices = itdmaxindex;
envindices = envmaxindex;
ildindices = ild;
