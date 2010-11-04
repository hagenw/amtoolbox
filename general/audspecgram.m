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
%               {'ihc','hilbert', 'noadapt', 'nomf'}
%
%-   'mlp',f   - Modulation low-pass filter to frequency f. Default is to
%                low-pass filter to 50 Hz.
%
%-   'nomf'    - No modulation filtering of any kind.
%
%-   'gammatonedelay' - Default
%
%-   'zerodelay' - New.
%
%-   'image'   - Use 'imagesc' to display the spectrogram. This is the default.
%
%-   'clim',[clow,chigh] - Use a colormap ranging from clow to chigh. These
%               values are passed to IMAGESC. See the help on IMAGESC.
%
%-   'dynrange',r - Limit the displayed dynamic range to r. This option
%                is especially usefull when displaying on a dB scale (no adaptation).
%
%-   'fullrange' - Use the full dynamic range. This is the default.
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
%-  'colorbar' - Display the colorbar. This is the default.
%
%-  'nocolorbar' - Do not display the colorbar.
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

% Define initial value for flags and key/value pairs.
definput.flags.adapt={'adapt','noadapt'};
definput.flags.thr={'nothr','thr'};
definput.flags.dynrange={'fullrange','dynrange'};
definput.flags.plottype={'image','contour','mesh','surf'};
definput.flags.clim={'noclim','clim'};
definput.flags.fmax={'nofmax','fmax'};
definput.flags.mlp={'mlp','nomf'};
definput.flags.delay={'gammatonedelay','zerodelay'};

definput.keyvals.ihc='dau';
definput.keyvals.dynrange=100;
definput.keyvals.thr=0;
definput.keyvals.clim=[0,1];
definput.keyvals.fmax=0;
definput.keyvals.ytick=[0,100,250,500,1000,2000,4000,8000];
definput.keyvals.mlp=50;
definput.keyvals.frange=[0,8000];
definput.keyvals.xres=800;
definput.keyvals.yres=600;

definput.flags.colorbar={'colorbar','nocolorbar'};

definput.groups.classic={'ihc','hilbert', 'noadapt', 'nomf'};


[flags,keyvals]=ltfatarghelper({},definput,varargin);

siglen=length(insig);

fhigh=keyvals.frange(2);
flow =keyvals.frange(1);

audlimits=freqtoaud('erb',keyvals.frange);

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
fc = erbspace(flow,fhigh,keyvals.yres);

% Calculate filter coefficients for the gammatone filter bank.
[gt_b, gt_a, delay]=gammatone(fc, fs, 'complex');

% Apply the Gammatone filterbank
outsig = 2*real(filterbank(gt_b,gt_a,insig,hopsize));

% The subband are now (possibly) sampled at a lower frequency than the
% original signal.
fssubband=round(fs/hopsize);

outsig = ihcenvelope(outsig,fssubband,keyvals.ihc);

if flags.do_adapt
  % non-linear adaptation loops
  outsig = adaptloop(outsig, fssubband);
end;
  
if flags.do_mlp
  % Calculate filter coefficients for the 50 Hz modulation lowpass filter.
  mlp_a = exp(-keyvals.mlp/fssubband);
  mlp_b = 1 - mlp_a;
  mlp_a = [1, -mlp_a];
  
  % Apply the low-pass modulation filter.
  outsig = filter(mlp_b,mlp_a,outsig);
end;
  
if flags.do_thr
  % keep only the largest coefficients.
  outsig=largestr(outsig,keyvals.thr);
end
  
% Apply transformation to coefficients.
if flags.do_noadapt
  % This is a safety measure to avoid log of negative numbers.
  outsig(:)=max(outsig(:),eps);

  outsig=20*log10(outsig);
end;

% 'dynrange' parameter is handled by threshholding the coefficients.
if flags.do_dynrange
  maxclim=max(outsig(:));
  outsig(outsig<maxclim-keyvals.dynrange)=maxclim-keyvals.dynrange;
end;

% Set the range for plotting
xsamples=siglen/hopsize;
xr=(0:hopsize:siglen-1)/fs;
yr=linspace(audlimits(1),audlimits(2),keyvals.yres);

% Determine the labels and position for the y-label.
ytickpos=freqtoerb(keyvals.ytick);

if flags.do_zerodelay
  % Correct the delays
  for n=1:keyvals.yres
    cut=round(delay(n)*fssubband);
    %size(outsig(:,n))
    %xsamples
    %cut
    %size([outsig(cut:end,n);zeros(cut-1,1)])
    outsig(:,n)=[outsig(cut:end,n);zeros(cut-1,1)];
  end;
end;

% Flip the output correctly. Each column is a subband signal, and should
% be display as the rows.
outsig=outsig.';

if flags.do_image
  if flags.do_clim
    imagesc(xr,yr,outsig,clim);
  else
    imagesc(xr,yr,outsig);
  end;
end;

if flags.do_contour
  contour(xr,yr,outsig);
end;

if flags.do_surf
  surf(xr,yr,outsig);
end;

if flags.do_mesh
  mesh(xr,yr,outsig);
end;

set(gca,'YTick',ytickpos);
% Use num2str here explicitly for Octave compatibility.
set(gca,'YTickLabel',num2str(keyvals.ytick(:)));

axis('xy');
xlabel('Time (s)')
ylabel('Frequency (Hz)')

if flags.do_colorbar
  colorbar;
end;

if nargout>0
  varargout={outsig,fc};
end;