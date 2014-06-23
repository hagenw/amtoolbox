function out = plotbaumgartner2014( p,tang,rang,varargin)
%PLOTBAUMGARTNER2014 plot probabilistic prediction matrixes
%   Usage:    plotbaumgartner2014(p,tang,rang);
%             plotbaumgartner2014(p,tang,rang,exptang,exprang);
%
%   Input parameters:
%     p       : prediction matrix containing probability mass vectors (PMVs) 
%               for the polar response angle as a function of the polar  
%               target angle (1st dim: response angle, 2nd dim: target
%               angle)
%     rang    : polar response angles
%     tang    : polar target angles
%
%   `plotbaumgartner2014(p,rang,tang)` plots predicted PMVs referring to  
%   the polar response angles *rang* as a function of the target angles 
%   *tang* with gray color coded probabilities similar to Baumgartner et al. 
%   (2002). Actual response patterns from psychoacoustic experiments can be 
%   overlayed optionally.
%
%   `h=plotbaumgartner2014(...)` additionally returns the figure handle.
%
%   `plotbaumgartner2014` accepts the following optional parameters:
%
%     'exptang',exptang   Overlay actual response patterns with the   
%                         experimetal polar target angles *exptang*.
%
%     'exprang',exprang   Experimetal polar response angles *exprang*
%                         corresponding to *exptang*.
%
%     'MarkerSize',ms     Set the marker (circles) size of the overlaying
%                         response pattern to *ms*. Default value is 6.
%
%     'cmax',cmax         Set the maximum probability of the color code to 
%                         *cmax*. Default value is 0.2.
%
%   `plotbaumgartner2014` takes the following flags at the end of the line 
%   of input arguments:
%
%     'colorbar'    Display the colorbar. This is the default.
%    
%     'nocolorbar'  Do not display the colorbar.
%
%   See also: baumgartner2014
%
%   References: baumgartner2014
  
% AUTHOR : Robert Baumgartner

definput.keyvals.exptang = [];
definput.keyvals.exprang = [];
definput.keyvals.MarkerSize = 6;
definput.keyvals.cmax = 0.1;
definput.flags.colorbar = {'colorbar','nocolorbar'};
[flags,kv]=ltfatarghelper({'exptang','exprang','MarkerSize','cmax'},definput,varargin);


%% Error handling
if size(p,1) ~= length(rang) || size(p,2) ~= length(tang)
  fprintf('\n Error: Dimension mismatch between p and rang/tang! \n Check the order of input arguments to fit plotbaumgartner2014(p,tang,rang). \n')
  return
end



%% Plot prediction matrix 
h = pcolor(tang,rang,p);
if not(isoctave) % does not work with x11 (mac osx)
  axis equal
end

set(gca,'XTick',-60:30:tang(end),'YTick',-60:30:rang(end),...    
	'XLim',[tang(1)-5,tang(end)+5],'YLim',[rang(1)-5,rang(end)+5],...
	'XMinorTick','on','YMinorTick','on')

colormap bone
shading flat
caxis([0 kv.cmax])
if flags.do_colorbar
    colorbar('southoutside')
end
xlabel('Target Angle (deg)')
ylabel('Response Angle (deg)')

%% Plot response pattern on top
if length(kv.exptang)==length(kv.exprang)
    hold on 
    h1 = plot( kv.exptang, kv.exprang, 'wo');  % shadow
    set(h1,'MarkerSize',kv.MarkerSize+round(kv.MarkerSize/3),'MarkerFaceColor','none') 
    h2 = plot( kv.exptang, kv.exprang, 'ko'); 
    set(h2,'MarkerSize',kv.MarkerSize,'MarkerFaceColor','none') 
    hold off
end

if nargout == 1
    out = h;
end
    
end