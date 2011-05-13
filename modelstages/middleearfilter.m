function b = middleearfilter(fs,varargin)
%MIDDLEEARFILTER   Middle ear filter.
%   Usage: b=middleearfilter(fs,varargin);
%          b=middleearfilter(fs);
%          b=middleearfilter;
%
%   MIDDLEEARFILTER(fs) computes the filter coefficients of a FIR
%   filter approximating the effect of the middle ear.
%
%   The following parameter and flags can be specified additionally:
%
%-    'order',order - Sets the filter order of the computed FIR filter.
%                     Default value is 512.
%
%-    'minimum' - calculates a minimum phase filter. This is the default.
%
%-    'zero' - returns a filter with zero phase. Since Matlab shifts the
%                symmetric impulse response due to no negative indices.
%                This results in a linear phase and hence a delay in the 
%                signal chain.
%
%   MIDDLEEARFILTER without any input arguments returns a table
%   describing the frequency response of the middle ear filter. First
%   column of the table contain frequencies and the second column
%   contains the amplitude (stapes peak velocity in m/s at 0dB SPL) of the 
%   frequency like in figure 2b) of Lopez-Poveda and Meddis 2001
%
%
%R  goode1994nkf lopezpoveda2001hnc
%
%   AUTHOR: Peter L. Soendergaard, Katharina Egger

%% ------ Check input options --------------------------------------------

% Define input flags
definput.flags.phase = {'minimum','zero'};
definput.keyvals.order = 512;

% Parse input options
[flags,keyvals]  = ltfatarghelper({},definput,varargin);

%%

data = data_lopezpoveda2001('fig2b', 'noplot');

if nargin==0
    b = data;
else
    if fs<=20000
        % In this case, we need to cut the table because the sampling
        % frequency is too low to accomodate the full range.
        indx=find(data(:,1)<fs/2);
        data = data(1:indx(end),:);
    else
        % otherwise the table will be extrapolated towards fs/2
        % data point added every 1000Hz
        lgth = size(data,1);
        for ii = 1:floor((fs/2-data(end,1))/1000)
            data(lgth+ii,1) = data(lgth+ii-1,1) + 1000;
            data(lgth+ii,2) = data(lgth+ii-1,2) / 1.1; % 1.1 corresponds to the decay of the last amplitude values
                                                  % = approx. ratio between amplitudes of frequency values seperated by 1000Hz                                 
        end
    end;  
    
    % for the function fir2 the last data point has to be at fs/2
    lgth = size(data,1);
    if data(lgth,1) ~= fs/2
        data(lgth+1,1) = fs/2;
        data(lgth+1,2) = data(lgth,2) / (1+(fs/2-data(lgth,1))*0.1/1000);
    end
    
    % Extract the frequencies and amplitudes, and put them in the format
    % that fir2 likes.
    freq=[0;...
        data(:,1).*(2/fs);...
       ];
    ampl=[0;...
        data(:,2);...
       ];

    b = fir2(keyvals.order,freq,ampl);
    
    b = b / 20e-6;      % scaling for SPL in dB re 20uPa
    
    if flags.do_minimum
        X = fft(b);
        Xmin = abs(X) .* exp(-i*imag(hilbert(log(abs(X)))));
        b = real(ifft(Xmin));        
    end
    
end

% if flags.do_plot
%     % Manually calculate the frequency response
%     fmid = abs(fftreal(b));
%     % Half the filter length.
%     n2=length(fmid);
%     % x-values for plotting.
%     xplot=linspace(0,fs/2,n2);
%     loglog(xplot/1000,fmid);
%     xlabel('Frequency (kHz)');
%     ylabel('FIR middleearfilter');
% end
