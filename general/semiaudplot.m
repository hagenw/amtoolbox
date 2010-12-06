function h=semiaudplot(x,y,varargin)
%SEMIAUDPLOT  2D plot on auditory scale
%   Usage: h=semiaudplot(x,y);
%
%   SEMIAUDPLOT(x,y) plots the data (x,y) on an auditory scale. By
%   default the values of the x-axis will be show on the Erb-scale.
%
%   SEMIAUDPLOT takes the following parameters at the end of the line of input
%   arguments:
%
%      'x'    - Make the x-axis use the auditory scale. This is the default.
%
%      'y'    - Make the y-axis use the auditory scale.
%
%      'erb'  - Use the Erb-scale. This is the default.
%
%      'bark' - Use the bark-scale.
%
%      'mel'  - Use the mel-scale.
%
%      'opts',c - Pass options stored in a cell array onto the plot function.
%    
%
%   See also: freqtoaud

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;


definput.flags.plotdir= {'x','y'};
definput.flags.scale  = {'erb','bark','mel'};
definput.keyvals.tick = [0,100,250,500,1000,2000,4000,8000];
definput.keyvals.res  = 500;
definput.keyvals.opts = {};

[flags,keyvals]=ltfatarghelper({},definput,varargin);

  
  
tick=[0,100,250,500,1000,2000,4000,8000];
n=500;
tickpos=freqtoaud(keyvals.tick,flags.scale);
  
  
if flags.do_x
  xmin=min(x);
  xmax=max(x);
  audminmax=freqtoaud([xmin,xmax],flags.scale);
  
  plotval=spline(x,y,audspace(xmin,xmax,n,flags.scale));
  plot(linspace(audminmax(1),audminmax(2),n),plotval,keyvals.opts{:});
  set(gca,'XTick',tickpos);
  set(gca,'XTickLabel',num2str(tick(:)));
 
end;

if flags.do_y
  ymin=min(y);
  ymax=max(y);
  audminmax=freqtoaud([ymin,ymax],flags.scale);
  
  plot(x,freqtoerb(y),keyvals.opts{:});
  set(gca,'YTick',tickpos);
  set(gca,'YTickLabel',num2str(tick(:)));
  
end;


