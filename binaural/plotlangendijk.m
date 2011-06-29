function [ out ] = plotlangendijk( p,rang,tang,name,latang,delta,deltat,cond,colbar )
% PLOTLANGENDIJK plots pdf-matrixes with gray colormap according to
% Langendijk et al. (2002)
% Usage:    plotlangendijk(p,rang,tang)
%           plotlangendijk(p,rang,tang,name)
%           plotlangendijk(p,rang,tang,name,latang,delta,deltat,cond)
%           plotlangendijk(p,rang,tang,name,latang,delta,deltat,cond,colbar)
% Input arguments:
%     p:        pdf-matrix for all target and response positions
%     rang:     response angles
%     tang:     target angles
%     name,latang,delta,deltat,cond: only for plot title
%     colbar:   switch for plotting colorbar-legend, use string 'colorbar' 
%               for switching on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR : Robert Baumgartner, OEAW Acoustical Research Institute
% latest update: 2010-08-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure('Name','Localization model','NumberTitle','off');
% clf
pcolor(tang,rang,p)
load('LangendijkColormap','langecmap');
set(gcf,'Colormap',langecmap);
shading interp
axis square

xlabel('Target Angle (°)')
ylabel('Response Angle (°)')
if ~exist('name','var')
    name='localization model';
end
if ~exist('latang','var')
    title(name)
else
    title(sprintf('condition: %s;  name: %s;  lateral angle: %d° +- %s°;  target: +- %s°',cond,name,latang,num2str(delta/2,'%7.1f'),num2str(deltat/2,'%7.1f')))
end

% plot legend
if ~exist('colbar','var')
    colbar='no';
end
if strcmp(colbar, 'colorbar')==1
   colorbar
end

out=gcf;
end

