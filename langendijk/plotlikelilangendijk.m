function plotlikelilangendijk(la,le,ci,tit)
% PLOTLIKELILANGENDIJK plots likelihood statistics according to Langendijk
% et al. (2002)
% Usage:           plotlikelilangendijk(la,le,ci)
%                  plotlikelilangendijk(la,le,ci,tit)
% Input arguments:
%     la:          actual likelihood
%     le:          expected likelihood
%     ci:          confidence interval for expected likelihood
%     [tit:]       set to 'spatstrat' optionally
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Robert Baumgartner, OEAW Acoustical Research Institute
% latest update: 2010-08-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('tit','var')
    tit='';
end
figure('Name','Likelihood statistic','NumberTitle','off');
clf
hold on
plot(0.5:length(la)+0.5,275*ones(length(la)+1,1),'k:') % unimodal
plot(0.5:length(la)+0.5,350*ones(length(la)+1,1),'k:') % bimodal
plot(0.5:length(la)+0.5,400*ones(length(la)+1,1),'k:') % trimodal
plot(0.5:length(la)+0.5,437*ones(length(la)+1,1),'k:') % flat
h=bar(la);
set(gca,'XTick',1:length(la))
if strcmp('spatstrat',tit)==1
    set(gca,'XTickLabel',{'baseline';'dummy';'warped'})
end
set(gca,'YLim',[200 550],'Layer','top')
set(gca,'Box','on')
set(h,'FaceColor','white','BarWidth',0.5)
ylabel('Likelihood')
xlabel('Condition')
errorbar(le,(ci(:,1)-ci(:,2))/2,'k.');
hold off

end