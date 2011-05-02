function data = data_lopezpoveda2001(varargin)
%DATA_LOPEZPOVEDA2001 Returns data points from the Lopez-Poveda and Meddis (2001) paper
%   Usage: data = data_lopezpoveda2001(flag)
%
%   DATA_LOPEZPOVEDA2001(flag) returns data points from the Lopez-Poveda and Meddis (2001)
%   The flag may be one of:
%
%-    'noplot'       - Don't plot, only return data. This is the default.
%
%-    'plot'         - Plot the data.
%  
%-    'fig2a'        - Data from Fig. 2(a), outer ear filter.
%
%-    'fig2b'        - Data from Fig. 2(b), middle ear filter.
%
%-    'fig2'         - Only for plotting: fig2a and fig2b combined
%
%R  lopezpoveda2001hnc

%   AUTHOR: Peter L. Soendergaard, Katharina Egger


%% ------ Check input options --------------------------------------------

% Define input flags
definput.flags.plot = {'noplot','plot'};
definput.flags.type = {'fig2a','fig2b'};

% Parse input options
[flags,keyvals]  = ltfatarghelper({},definput,varargin);

%% ------ Data points from the paper ------------------------------------
%
% Data for the given figure

if flags.do_fig2a

  data=data_pralong1996;
  
  if flags.do_plot

    fs=22050;
    bout=headphonefilter(fs);

    % Manually calculate the frequency response.
    fout = 20*log10(abs(fftreal(bout)));

    % Half the filter length.
    n2=length(fout);

    figure;
    hold on;
    % Plot the measured data
    x=data(:,1);
    freqresp=20*log10(data(:,2));
    semilogx(x,freqresp,'ro');
    
    % Plot the filter
    x_filter=linspace(0,fs/2,n2);
    semilogx(x_filter,fout);
    hold off;
  end;
  
end;

if flags.do_fig2b
  % get peak-to-peak displacement data
  data = data_goode1994;

  % get velocity (proportional voltage) acc. to formula in Goode et al. (1994), page 147
  data(:,2) = data(:,2) * 1e-6 * 2 * pi .* data(:,1);  
  
  % to get data at 0dB SPL (assumed that stapes velocity is linearly related to pressure
  data(:,2) = data(:,2) * 10^(-104/20);
  
  % to get stapes PEAK velocity, multiply amplitudes by sqrt(2)
  data(:,2) = data(:,2).*sqrt(2);   
  
  
  % extrapolated data points, directly read from figure 2b) of Lopez-Poveda
  % and Meddis (2001)
  extrp = [100 1.181E-09; ...
           200	2.363E-09; ...
           7000 8.705E-10; ...
           7500 8.000E-10; ...
           8000 7.577E-10; ...
           8500 7.168E-10; ...
           9000 6.781E-10; ...
           9500 6.240E-10; ...
           10000 6.000E-10];
  
  if flags.do_plot
    figure;
    loglog(data(:,1)/1000,data(:,2),'ok', 'MarkerFaceColor', 'k');
    hold on
    loglog(extrp(:,1)/1000,extrp(:,2),'ok');
    xlabel('Frequency (kHz)');
    ylabel('Stapes peak velocity (m/s) at 0dB SPL)');
    axis([0.1,10,1e-10,1e-7]);
  end;

  data = [extrp(extrp(:,1) < data(1,1),:); data; extrp(extrp(:,1) > data(end,1),:)];

end;

