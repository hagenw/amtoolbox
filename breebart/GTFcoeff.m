function [a,b] = GTFcoeff(in,style,fs)

% Calculate the filter coefficients for the gammatone filterbank
%
% Usage: [a,b] = GTFcoeff(in,style,fs)
%
%   a - feedback coefficients
%   b - feedforward coefficients
%
%   in - array of filter center frequencies (can be single value)
%   style - can be '3rdorder' or '4thorder'
%   fs - array of sampling rates corresponding to the filter center
%       frequencies
%
% rev. 0.01 2006-Oct-05

switch style
    case '3rdOrder' % based on Hohmann (2002), impulse invariant all-pole design
        order = 3;
        alpha = pi * factorial(2 * order - 2) * 2^(-(2 * order - 2)) / factorial(order - 1)^2; %Hohmann, 2002
        beta = 1/alpha * (24.7 + in/9.265);  
        b = zeros(1,length(in));
        a = zeros(order+1,length(in));
        for idx = 1:length(in)
            [bt,at] = gammaf1(in(idx),beta(idx),order,2,fs(idx));
            b(:,idx)=conj(bt');
            a(:,idx)=conj(at');
        end

    case '4thOrder' % based on Hohmann (2002), impulse invariant all-pole design
        order = 4;
        alpha = pi * factorial(2 * order - 2) * 2^(-(2 * order - 2)) / factorial(order - 1)^2; %Hohmann, 2002
        beta=1/alpha * (24.7 + in/9.265);
        b = zeros(1,length(in));
        a = zeros(order+1,length(in));
        for idx=1:length(in)
            [bt,at]=gammaf1(in(idx),beta(idx),4,2,fs(idx));
            b(:,idx)=conj(bt');
            a(:,idx)=conj(at');
        end

    otherwise
        error(['GTFcoeff: style ''' style ''' unknown'])
end


% erb = 0.14*CF;
% beta = 1.0184*erb;      %0.637 for 2nd order
% [b,a] = gammaf1(CF,beta,4,2,fs);
% out = 2*real(filter(b,a,in));

%eof