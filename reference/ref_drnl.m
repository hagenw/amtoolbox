function out = ref_drnl(x,CF,fs)
% script for pemo preprocessing using an implementation of the dual
% resonance nonlinear (DRNL) filter (Lopez-Poveda, meddis 2001)
% The filter models the BM non-linearity
% Author: Morten Løve Jepsen, 2.nov 2005
%
% usage: out = drnl(x,CF,fs)

[linDRNLpar,nlinDRNLpar] = getDRNLparam(CF);

[GTlin_b,GTlin_a] = coefGtDRNL(linDRNLpar(1).vals,linDRNLpar(3).vals,1,fs); %get GT filter coeffs
[LPlin_b,LPlin_a] = coefLPDRNL(linDRNLpar(5).vals,fs); % get LP filter coeffs



%% Apply the middle-ear filter
me_fir = middleearfilter(fs);
x = filter(me_fir,1,x);

% linDRNLpar(4).vals
% pause
y_lin = x.*linDRNLpar(4).vals; % Apply linear gain

% Now filtering
for n = 1:linDRNLpar(2).vals % Gammatone filtering multiple times for cascading
    y_lin = real(filter(GTlin_b,GTlin_a,y_lin));  
end
for n = 1:linDRNLpar(6).vals % cascade of lowpass filters
    y_lin = filter(LPlin_b,LPlin_a,y_lin);           
end
% end of linear part %%%%%%%%%%%%%%%%%%%%%%%

% Non-linear part%%%%%%%%%%%%%%%%%%%%%%%%%%%
[GTnlin_b,GTnlin_a] = coefGtDRNL(nlinDRNLpar(1).vals,nlinDRNLpar(3).vals,1,fs); %get GT filter coeffs
[LPnlin_b,LPnlin_a] = coefLPDRNL(nlinDRNLpar(7).vals,fs); % get LP filter coeffs

y_nlin = x;

% Now GT filtering
for n = 1:nlinDRNLpar(2).vals % Gammatone filtering multiple times for cascading
    y_nlin = real(filter(GTnlin_b,GTnlin_a,y_nlin));  
end

% Broken stick nonlinearity
a = nlinDRNLpar(4).vals;
b = nlinDRNLpar(5).vals;
c = nlinDRNLpar(6).vals;
y_decide = [a*abs(y_nlin); b*(abs(y_nlin)).^c];
y_nlin = sign(y_nlin).* min(y_decide);

% Now GT filtering again
for n = 1:nlinDRNLpar(2).vals % Gammatone filtering multiple times for cascading
    y_nlin = real(filter(GTnlin_b,GTnlin_a,y_nlin));  
end
% then LP filtering
if nlinDRNLpar(8).vals > 0
for n = 1:nlinDRNLpar(8).vals % cascade of lowpass filters
    y_nlin = filter(LPnlin_b,LPnlin_a,y_nlin);           
end
end

out = (y_lin + y_nlin);
