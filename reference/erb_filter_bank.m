function y = erb_filter_bank(x, filter_coefficients)
%__________________________________________________________________________
%
% y = erb_filter_bank(x, filter_coefficients)
%
% Process an input waveform with a gammatone filter bank. This function 
% takes a single sound vector, and returns an array of filter outputs, one 
% channel per row.
%
% The filter_coefficients parameter, which completely specifies the 
% Gammatone filterbank, should be designed with the 
% make_erb_filter_coefficients function.  If it is omitted, the filter 
% coefficients are computed for you assuming a 22050Hz sampling rate and 
% 64 filters regularly spaced on an ERB scale from fs/2 down to 100Hz.
%
% Options:
%   x                       - waveform signal
%   filter_coefficients     - filter coefficients for the erb filter bank
%
% Outputs:
%   y                       - the filtered waveform signal
%__________________________________________________________________________
%bi.mo                                                                v.0.3

%
% Malcolm Slaney @ Interval, June 11, 1998.
% (c) 1998 Interval Research Corporation  
% Thanks to Alain de Cheveigne' for his suggestions and improvements.
%
% changes by:
% Hagen Wierstorf
% T-Labs, Berlin
% hagen.wierstorf@telekom.de
% 2009/05/06
%


%% Make global variables/constants avaiable
global def


%% Check arguments
if nargin < 1
	error('Syntax: output_array = erb_filter_bank(input_vector[, filter_coefficients]);');
end

if nargin < 2
	filter_coefficients = make_erb_filter_coefficients(def.fs,64,100);
end

if size(filter_coefficients,2) ~= 10
	error('filter_coefficients parameter passed to erb_filter_bank is the wrong size.');
end

if size(x,2) < size(x,1)
	x = x';
end


%% Extract/name the given filter coefficients
A0  = filter_coefficients(:,1);
A11 = filter_coefficients(:,2);
A12 = filter_coefficients(:,3);
A13 = filter_coefficients(:,4);
A14 = filter_coefficients(:,5);
A2  = filter_coefficients(:,6);
B0  = filter_coefficients(:,7);
B1  = filter_coefficients(:,8);
B2  = filter_coefficients(:,9);
gain= filter_coefficients(:,10);	


%% Filter and generate output signal
output = zeros(size(gain,1), length(x));
for chan = 1:size(gain,1)
	% What does this filter cycle do?
    y1 = filter([A0(chan)/gain(chan) A11(chan)/gain(chan) ...
		         A2(chan)/gain(chan)], ...
				[B0(chan) B1(chan) B2(chan)], x);
	y2 = filter([A0(chan) A12(chan) A2(chan)], ...
				[B0(chan) B1(chan) B2(chan)], y1);
	y3 = filter([A0(chan) A13(chan) A2(chan)], ...
				[B0(chan) B1(chan) B2(chan)], y2);
	y4 = filter([A0(chan) A14(chan) A2(chan)], ...
				[B0(chan) B1(chan) B2(chan)], y3);
    y(chan,:) = y4;
end


%% Debug
%if 0
%	semilogx((0:(length(x)-1))*(fs/length(x)),20*log10(abs(fft(output))));
%end


%% EOF
