function varargout=audspecgram(insig,fs,varargin)
%AUDSPECGRAM  Auditory spectrogram.
%   Usage: audspecgram(insig,fs,op1,op2, ... );
%          C=audspecgram(insig,fs, ... );
%
%   AUDSPECGRAM(insig,fs) plots an auditory spectrogram of the signal insig,
%   which has been sampled at a sampling rate of fs Hz. The output is
%   low-pass modulation filtered before presentation.
%
%   The frequency axis is diplayed on a erb-scale, but labelled in
%   Hz. Using the mouse to get plot coordinates will reveal the real
%   value in erb's. Use ERBTOFREQ to convert to Hz.
%
%   C=AUDSPECGRAM(insig,fs, ... ) returns the image to be displayed as a
%   matrix. Use this in conjunction with IMWRITE etc. Do NOT use this as a
%   method to compute an auditory representation. Use some of the model
%   preprocessing functions for this.
%
%   Be carefull with long signals, as the routine may lock up the
%   interpreter.
%
%   Additional arguments can be supplied like this:
%   AUDSPECGRAM(insig,fs,'dynrange',30). The arguments must be character
%   strings possibly followed by an argument:
%
%-   'adapt'   - Model adaptation. This is the default. This options also
%                sets the output to be displayed on a linear scale.
%
%-   'noadapt' - Do not model adaptation. This option also sets a Db scale to
%                display the output.
%
%-   'ihc',modelname - Pass modelname to IHCENVELOPE to determine the inner
%                hair cell envelope extraction process to use. Default is to
%                use the 'dau' model.
%
%-   'classic' - Display a classic spectrogram. This option is equal to
%               'hilbert', 'noadapt', 'nomf'
%
%-   'mlp',f   - Modulation low-pass filter to frequency f. Default is to
%                low-pass filter to 50 Hz.
%
%-   'nomf'    - No modulation filtering of any kind.
%
%-   'image'   - Use 'imagesc' to display the spectrogram. This is the default.
%
%-   'clim',[clow,chigh] - Use a colormap ranging from clow to chigh. These
%               values are passed to IMAGESC. See the help on IMAGESC.
%
%-   'dynrange',r - Limit the displayed dynamic range to r. This option
%                is especially usefull when displaying on a log/Db scale.
%
%-   'fullrange' - Use the full dynamic range.
%
%-   'ytick'   - A vector containing the frequency in Hz of the yticks.
%
%-   'thr',r   - Keep only the largest fraction r of the coefficients, and
%                set the rest to zero.
%
%-   'frange',[flow,fhigh] - Choose a frequency scale ranging from flow to
%                fhigh, values are eneterd in Hz. Default is to display from
%                0 to 8000 Hz.
%
%-   'xres',xres - Approximate number of pixels along x-axis / time.
%
%-   'yres',yres - Approximate number of pixels along y-axis / frequency If
%                 only one of 'xres' and 'yres' is specified, the default
%                 aspect ratio will be used.
%
%-   'displayratio',r - Set the default aspect ratio.
%
%-   'contour' - Do a contour plot to display the spectrogram.
%          
%-   'surf'    - Do a surf plot to display the spectrogram.
%
%-   'mesh'    - Do a mesh plot to display the spectrogram.
%
%-
%
%   See also:  erbtofreq, dau96
%
%   Demos:  demo_audspecgram
  
%   AUTHOR : Peter Soendergaard.
%   TESTING: NA
%   REFERENCE: NA
  
if nargin<2
  error('Too few input arguments.');
end;

if ~isnumeric(insig) || ~isvector(insig)
  error('%s: Input must be a vector.',upper(mfilename));
end;

global AMT_CONF;

if isempty(AMT_CONF)
  error(['You need to run AMTSTART. This could be because you have used a ', ...
        'CLEAR ALL statement.']);
end;

% Get the default values.
defaults=AMT_CONF.plotdefaults;

% Approximate resolution along time-axis.
xres=TF_CONF.xres;

% Ratio of frequency resolution versus time resolution
displayratio=TF_CONF.displayratio;

yres=ceil(xres*displayratio);

% Define initial value for flags and key/value pairs.
defnopos.flags.adapt={'adapt','noadapt'};
defnopos.flags.thr={'nothr','thr'};
defnopos.flags.dynrange={'fullrange','dynrange'};
defnopos.flags.plottype={'image','contour','mesh'};
defnopos.flags.clim={'noclim','clim'};
defnopos.flags.fmax={'nofmax','fmax'};
defnopos.flags.mlp={'mlp','nomf'};

defnopos.keyvals.ihc='dau';
defnopos.keyvals.dynrange=100;
defnopos.keyvals.thr=0;
defnopos.keyvals.clim=[0,1];
defnopos.keyvals.fmax=0;
defnopos.keyvals.ytick=[0,50,100,200,500,1000,2000,4000,8000];
defnopos.keyvals.mlp=50;


[flags,keyvals,fs]=ltfatarghelper(1,{[]},defnopos,varargin,upper(mfilename));

if ~doxres
  xres=floor(yres/displayratio);
end;

if ~doyres
  yres=ceil(xres*displayratio);
end;

siglen=length(insig);

fhigh=frange(2);
flow =frange(1);

audlimits=freqtoaud('erb',frange);

% fhigh can at most be the Nyquest frequency
fhigh=min(fhigh,fs/2);

% Downsample this signal if it is sampled at a much higher rate than
% 2*fhigh. This reduces memory consumption etc. 1.5 and 1.2 are choosen as a
% safeguard to not loose information.
if fs>2*1.5*fhigh
  
  fsnew=round(fhigh*2*1.2);

  % Determine new signal length
  siglen=round(siglen/fs*fsnew);
  
  % Do the resampling using an FFT based method, as this is more flexible
  % than the 'resample' method included in Matlab
  insig=fftresample(insig,siglen);

  % Switch to new value
  fs=fsnew;  
end;

% Determine the hopsize
% Using a hopsize different from 1 is currently not possible because all
% the subsequent filters fail because of a to low subband sampling rate.
%hopsize=max(1,floor(siglen/xres));

hopsize=1;

% find the center frequencies used in the filterbank
fc = erbspace(flow,fhigh,yres);

% Calculate filter coefficients for the gammatone filter bank.
[gt_b, gt_a]=gammatone(fc, fs);

% Apply the Gammatone filterbank
outsig = filterbank(gt_b,gt_a,insig,hopsize);

% The subband are now (possibly) sampled at a lower frequency than the
% original signal.
fssubband=round(fs/hopsize);

if dorectify && (fssubband>2000)
  outsig = ihcenvelope(2*real(outsig),fssubband,'dau');
else
  outsig = abs(outsig);
end;

if doadapt
  % non-linear adaptation loops
  outsig = adaptloop(outsig, fssubband);
end;
  
if domlp
  % Calculate filter coefficients for the 50 Hz modulation lowpass filter.
  mlp_a = exp(-mlp/fssubband);
  mlp_b = 1 - mlp_a;
  mlp_a = [1, -mlp_a];
  
  % Apply the low-pass modulation filter.
  outsig = filter(mlp_b,mlp_a,outsig);
end;
  
if domask
  % keep only the largest coefficients.
  outsig=largestr(outsig,mask_val);
end
  
% Apply transformation to coefficients.
if dolog
  % This is a safety measure to avoid log of negative numbers if the
  % users chooses an incorrect combination of options (e.g. 'adapt' and 'log').
  outsig(:)=max(outsig(:),realmin);

  outsig=20*log10(outsig);
end;

% 'dynrange' parameter is handled by threshholding the coefficients.
if dorange
  maxclim=max(outsig(:));
  outsig(outsig<maxclim-range)=maxclim-range;
end;

% Set the range for plotting
xr=(0:hopsize:siglen-1)/fs;
yr=linspace(audlimits(1),audlimits(2),length(fc));

% Determine the labels and position for the y-label.
ytickpos=freqtoerb(ytick);

% Flip the output correctly. Each column is a subband signal, and should
% be display as the rows.
outsig=outsig.';

switch(plottype)
  case 'image'
    if doclim
      imagesc(xr,yr,outsig,clim);
    else
      imagesc(xr,yr,outsig);
    end;
  case 'contour'
    contour(xr,yr,outsig);
  case 'surf'
    surf(xr,yr,outsig);
end;
set(gca,'YTick',ytickpos);
% Use num2str here explicitly for Octave compatibility.
set(gca,'YTickLabel',num2str(ytick(:)));

axis('xy');
xlabel('Time (s)')
ylabel('Frequency (Hz)')

if nargout>0
  varargout={outsig,fc};
end;