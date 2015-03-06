function [c,lags]=ncorr(x,y,numlags);


% function [c,lags]=ncorr(x,y,numlags);
%
% normalized cross-correlation function
%
% Jonas Braasch
% (c) 2013, Rensselaer Polytechnic Institute
% contact: jonasbraasch@gmail.com
%
% INPUT PARAMETERS:
% x       = Signal 1
% y       = Signal 2
% numlags = number of lags (tau) that will be analyzed
% OUTPUT PARAMETERS:
% c       = normalized cross-correlation
% lags    = lag coefficients 

[c,lags]=xcorr(x,y,numlags);
c=c./(sqrt(sum(x.^2)*sum(y.^2)));



