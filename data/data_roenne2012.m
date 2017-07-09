function [ur,fs]  = data_roenne2012(varargin)
%DATA_ROENNE2012 Unitary response
%   Usage: ur = data_roenne2012()
%
%   `[ur,fs]=data_roenne2012` returns the unitary response from Roenne
%   (2012) and its sampling frequency, $fs=30000$.
%
%   Examples:
%   ---------
%
%   The first plot shows the unitary response in the time-domain:::
%
%     [ur,fs]  = data_roenne2012;
%     plot((0:length(ur)-1)/fs,ur);
%     xlabel('Time / seconds');
%     ylabel('Amplitude');
%
%   The second plot shows the magnitude response of the unitary
%   response, normalized so the highest peak reaches 0-dB:::
%
%     [ur,fs]  = data_roenne2012;
%     magresp(ur,fs,90,'fir','1');
%
%   References: roenne2012modeling

% TODO: explain "unitary response"
  
s=amt_load('roenne2012','ur.mat');
ur=s.ur;

fs=30000;


