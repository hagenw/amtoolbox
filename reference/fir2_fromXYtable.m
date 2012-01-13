% function [b] = fir2_fromXYtable (xytable, fs, dBorAmp)
%
% Calculates the n+1 coefficients for an nth-order linear FIR filter from the
% frequency response specified in 'xytable' and a sampling freq. 'fs'.
% xytable is assumed to be such that the frequencies are specified in the
% first row (xytable(:,1)) and the magnitudes in the second (xytable(:,2)).
%		Parameter 'dBorAmp' specifies whether the magnitude data in the
% XY table must be interpreted as dB (gain) or as amplitude.
% The option 'dB' calculates amplitude values as:
% y = 10^(x/20).
% ------------------------------------------------------------------
function [b] = fir2_fromXYtable (order, xytable, fs, dBorAmp)

if max(xytable(:,1)) < (fs/2)
   [n,m] = size(xytable);
   
   % Assign frequencies.
   f_mag = zeros(n+2,m); 
   f_mag(1,1) = 0.0; % The first normalised freq. must be 0.0;
   f_mag(2:n+1,1) = xytable(:,1) .* (2.0 / fs);
   f_mag(n+2,1) = 1.0; %The last normalised freq. must be the Nyquist freq.
   
   % Assign amplitudes.
   f_mag(1,2) = 0.0; % The first normalised amplitude is set to a small value;
   switch dBorAmp
   case 'dB'
      f_mag(2:n+1,2) = dB2amp2(xytable(:,2));
   case 'amp'
      f_mag(2:n+1,2) = xytable(:,2);
   otherwise
      error('Invalid value for argument: dBorAmp.');
   end
   
   % Assign last amplitude value corresponding to Nyquist freq.
   f_mag(n+2,2) = 0.0;
   
   % Calculate FIR coefficients.
   b = fir2(order, f_mag(:,1), f_mag(:,2));
else
   error('Freq range in XY-table exceeds Nyquist freq (fs/2)');
end


%OLDFORMAT
