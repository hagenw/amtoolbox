% taken from PEMO framework, modified by MLJ 24. april 2007
%
% Usage: [inf, cerbs] = getGFBCenterERBs(lowf, uppf, density);
%
% in = input vector 
% lowf = lower frequency boundary for filters
% uppf = upper frequency boundary for filters
%        (if lowf == uppf, only one filter is chosen)
% fs = sampling rate in Hz
%
% out = output array (columns = channels, rows = time)
% cfs = vector of center frequencies in ERB (use erbtofreq.m to convert to Hz)

function [info, centerFrq] = getFBCenterERBs(lowf, uppf, baseF);

%baseF = 4000;

tmp = freqtoerb([lowf baseF uppf]);

if ( tmp(1) == tmp(3) )

	cerbs = tmp(1);

else 
	
	tmp2 = [0:100];
	tmp2 = [-fliplr(tmp2) tmp2(2:end)]+tmp(2);
	
	i_start = min(find(tmp2>=tmp(1)));
	i_end = max(find(tmp2<=tmp(3)));

	cerbs = tmp2(i_start:i_end);
end

info = length(cerbs);
centerFrq = erbtofreq(cerbs);

% eof

%OLDFORMAT
