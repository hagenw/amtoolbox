function filter_coefficients = make_erb_filter_coefficients( number_of_channels,lowest_frequency )
%__________________________________________________________________________
%
% [filter_coefficients] = make_erb_filter_coefficients( number_of_channels,lowest_frequency )
%
% This function computes the filter coefficients for a bank of 
% Gammatone filters. These filters were defined by Patterson and 
% Holdworth for simulating the cochlea.  
% 
% The result is returned as an array of filter coefficients.  Each row 
% of the filter arrays contains the coefficients for four second order 
% filters.  The transfer function for these four filters share the same
% denominator (poles) but have different numerators (zeros).  All of these
% coefficients are assembled into one vector that the function 
% erb_filterbank can take apart to implement the filter.
%
% The filter bank contains "number_of_channels" channels that extend from
% half the sampling rate (def.fs) to "lowest_frequency".  Alternatively, 
% if the number_of_channels input argument is a vector, then the values of 
% this vector are taken to be the center frequency of each desired filter.
% (The lowest_frequency argument is ignored in this case.)
%__________________________________________________________________________
%bi.mo                                                                v.0.3

%
% Note this implementation fixes a problem in the original code by
% computing four separate second order filters.  This avoids a big
% problem with round off errors in cases of very small center_frequencies 
% (100Hz) and large sample rates (44kHz).  The problem is caused by 
% roundoff error when a number of poles are combined, all very close to the 
% unit circle.  Small errors in the eigth order coefficient, are multiplied
% when the eigth root is taken to give the pole location.  These small
% errors lead to poles outside the unit circle and instability.  Thanks
% to Julius Smith for leading me to the proper explanation.
%
% Execute the following code to evaluate the frequency
% response of a 10 channel filterbank.
%	filter_coefficients = make_erb_filter_coefficients(16000,10,100);
%	y = ERBFilterBank([1 zeros(1,511)], filter_coefficients);
%	resp = 20*log10(abs(fft(y')));
%	freqScale = (0:511)/512*16000;
%	semilogx(freqScale(1:255),resp(1:255,:));
%	axis([100 16000 -60 0])
%	xlabel('Frequency (Hz)'); ylabel('Filter Response (dB)');
%
% Rewritten by Malcolm Slaney@Interval.  June 11, 1998.
% (c) 1998 Interval Research Corporation  
%
% changes by:
% Hagen Wierstorf
% T-Labs, Berlin
% hagen.wierstorf@telekom.de
% 2009/05/18
%

%% Make global variables/constants avaiable
global def


%% Check/calculate the center frequencies
% If only the number of channels are given, calculate the center
% frequencies,
if length(number_of_channels) == 1
	
    center_frequencies = erb_space(lowest_frequency, def.fs/2, number_of_channels);
    
% else use the given center frequencies.
else
    
	center_frequencies = number_of_channels(1:end);
    
    % Check the correct oriantation of the center frequencies array
	if size(center_frequencies,2) > size(center_frequencies,1)
		center_frequencies = center_frequencies';
    end
    
end


%% Calculate the ERB filter coefficients
T = 1/def.fs;
ERB = ((center_frequencies/def.erb_bank_ear_q).^def.erb_bank_order + def.erb_bank_minimum_bandwidth^def.erb_bank_order).^(1/def.erb_bank_order);
B = 1.019*2*pi*ERB;

A0 = T;
A2 = 0;
B0 = 1;
B1 = -2*cos(2*center_frequencies*pi*T)./exp(B*T);
B2 = exp(-2*B*T);

A11 = -(2*T*cos(2*center_frequencies*pi*T)./exp(B*T) + 2*sqrt(3+2^1.5)*T*sin(2*center_frequencies*pi*T)./ ...
		exp(B*T))/2;
A12 = -(2*T*cos(2*center_frequencies*pi*T)./exp(B*T) - 2*sqrt(3+2^1.5)*T*sin(2*center_frequencies*pi*T)./ ...
		exp(B*T))/2;
A13 = -(2*T*cos(2*center_frequencies*pi*T)./exp(B*T) + 2*sqrt(3-2^1.5)*T*sin(2*center_frequencies*pi*T)./ ...
		exp(B*T))/2;
A14 = -(2*T*cos(2*center_frequencies*pi*T)./exp(B*T) - 2*sqrt(3-2^1.5)*T*sin(2*center_frequencies*pi*T)./ ...
		exp(B*T))/2;

gain = abs((-2*exp(4*1i*center_frequencies*pi*T)*T + ...
             2*exp(-(B*T) + 2*1i*center_frequencies*pi*T).*T.* ...
             (cos(2*center_frequencies*pi*T) - sqrt(3 - 2^(3/2))* ...
               sin(2*center_frequencies*pi*T))) .* ...
           (-2*exp(4*1i*center_frequencies*pi*T)*T + ...
             2*exp(-(B*T) + 2*1i*center_frequencies*pi*T).*T.* ...
             (cos(2*center_frequencies*pi*T) + sqrt(3 - 2^(3/2)) * ...
               sin(2*center_frequencies*pi*T))).* ...
           (-2*exp(4*1i*center_frequencies*pi*T)*T + ...
             2*exp(-(B*T) + 2*1i*center_frequencies*pi*T).*T.* ...
            (cos(2*center_frequencies*pi*T) - ...
              sqrt(3 + 2^(3/2))*sin(2*center_frequencies*pi*T))) .* ...
           (-2*exp(4*1i*center_frequencies*pi*T)*T + 2*exp(-(B*T) + 2*1i*center_frequencies*pi*T).*T.* ...
             (cos(2*center_frequencies*pi*T) + sqrt(3 + 2^(3/2))*sin(2*center_frequencies*pi*T))) ./ ...
           (-2 ./ exp(2*B*T) - 2*exp(4*1i*center_frequencies*pi*T) +  ...
             2*(1 + exp(4*1i*center_frequencies*pi*T))./exp(B*T)).^4);
	
allfilts = ones(length(center_frequencies),1);

% filter coefficients
filter_coefficients = [A0*allfilts A11 A12 A13 A14 A2*allfilts B0*allfilts B1 B2 gain];


%% Debugging code
% if (0)						% Test Code
% 	A0  = filter_coefficients(:,1);
% 	A11 = filter_coefficients(:,2);
% 	A12 = filter_coefficients(:,3);
% 	A13 = filter_coefficients(:,4);
% 	A14 = filter_coefficients(:,5);
% 	A2  = filter_coefficients(:,6);
% 	B0  = filter_coefficients(:,7);
% 	B1  = filter_coefficients(:,8);
% 	B2  = filter_coefficients(:,9);
% 	gain= filter_coefficients(:,10);	
% 	chan=1;
% 	x = [1 zeros(1, 511)];
% 	y1=filter([A0(chan)/gain(chan) A11(chan)/gain(chan) ...
% 		A2(chan)/gain(chan)],[B0(chan) B1(chan) B2(chan)], x);
% 	y2=filter([A0(chan) A12(chan) A2(chan)], ...
% 			[B0(chan) B1(chan) B2(chan)], y1);
% 	y3=filter([A0(chan) A13(chan) A2(chan)], ...
% 			[B0(chan) B1(chan) B2(chan)], y2);
% 	y4=filter([A0(chan) A14(chan) A2(chan)], ...
% 			[B0(chan) B1(chan) B2(chan)], y3);
% 	semilogx((0:(length(x)-1))*(def.fs/length(x)),20*log10(abs(fft(y4))));
% end


