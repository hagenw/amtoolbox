function h=semiaudplot(plotdir,audscale,x,y,varargin)
%SEMIAUDPLOT  2D plot on auditory scale
%   Usage: h=semiaudplot(plotdir,audscale,x,y,varargin);
%
%   SEMIAUDPLOT(plotdir,audscale,x,y,varargin) plots the data (x,y) on an
%   auditory scale specified by audscale. If plotdir=='x' then the
%   x-scale is displayed on an auditory scale, and similarly with
%   plotdir=='y'.
%M
%   See also: freqtoaud
  
tick=[0,100,250,500,1000,2000,4000,8000];
n=500;
tickpos=freqtoaud(audscale,tick);
  
  
switch(lower(plotdir))
 case 'x'
  xmin=min(x);
  xmax=max(x);
  audminmax=freqtoaud(audscale,[xmin,xmax]);
  
  plotval=spline(x,y,erbspace(xmin,xmax,n));
  plot(linspace(audminmax(1),audminmax(2),n),plotval,varargin{:});
  set(gca,'XTick',tickpos);
  set(gca,'XTickLabel',num2str(tick(:)));
 
 case 'y'
  ymin=min(y);
  ymax=max(y);
  audminmax=freqtoaud(audscale,[ymin,ymax]);
  
  plot(x,freqtoerb(y),varargin{:});
  set(gca,'YTick',tickpos);
  set(gca,'YTickLabel',num2str(tick(:)));
  
 otherwise
  error('plotdir must be either x or y.');
end;

