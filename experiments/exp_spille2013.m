function example_output = exp_spille2013
%EXP_SPILLE2013 Experiments from Spille et al 2013. Model: Dietz et al 2011
%
%   `exp_spille2013` reproduce Fig. no. 6 from the Spille et
%   al. 2013 AABBA chapter.
%
%   Requires exp_spille2013_s123456.wav in amtbasepath/signals/
%
%   Examples:
%   ---------
%
%   To display Fig. 6 use :::
%
%     exp_spille2013;
%
%   See also: exp_dietz2011


% AUTHOR: Mathias Dietz

% choose input signal (e.g. spille2013_s123456)
fn = 'exp_spille2013_s123456.wav';
targetfn = fullfile(amtbasepath,'signals',fn);
if not(exist(targetfn,'file'))
  disp(['Downloading ' fn ' from http://amtoolbox.sourceforge.net/data/signals/']);  
  sourcefn = ['http://amtoolbox.sourceforge.net/data/signals/' fn];
  urlwrite(sourcefn,targetfn);
end
signal=wavread(targetfn);
fs=44100;
s_pos =[-75 -40 -10 10 40 75];

ic_threshold=0.98;
cn = [10 1]; % channel numbers for separate plots (1st entry also for time plot)
panellabel = 'acbd';

% run IPD model on signal
[hairc_fine, hairc_mod, fc, hairc_ild]=dietz2011(signal,fs);

% convert interaural information into azimuth
itd_unwrapped = ...
    dietz2011unwrapitd(hairc_fine.itd_lp,hairc_ild(:,1:12),hairc_fine.f_inst,2.5);
lookup = load('dietz2011itd2anglelookup.mat');
angl=itd2angle(itd_unwrapped,lookup);

h_ic=zeros(91,12);
h_all=histc(angl,-90:2:90);
for n=1:12
    h_ic(:,n)=histc(angl(hairc_fine.ic(:,n)>ic_threshold&[diff(hairc_fine.ic(:,n))>0; 0],n),-90:2:90);
end

% plot
figure;
fontsize = 14;
set(gcf,'Position',[100 100 1170 700])

for panel = 1:4
    subplot(2,2,panel)
    switch panel
        case 1
            bar(-90:2:90,h_all(:,cn(1)))
            hold on
            bar(-90:2:90,h_ic(:,cn(1)),'r')
            title(['Azimuth histogram of ch. ' num2str(cn(1)) ' at cf = ' ...
                num2str(round(fc(cn(1)))) ' Hz'],'Fontsize',fontsize)
            ymax = max(h_all(:,cn(1)));
        case 2
            bar(-90:2:90,h_all(:,cn(2)))
            hold on
            bar(-90:2:90,h_ic(:,cn(2)),'r')
            title(['Azimuth histogram of ch. ' num2str(cn(2)) ' at cf = ' ...
                num2str(round(fc(cn(2)))) ' Hz'],'Fontsize',fontsize)
            ymax = max(h_all(:,cn(2)));
        case 3
            t=(1:length(signal))*1/fs;
            plot(t,angl(:,cn(1)),'b.');
            hold on;
            plot(t(hairc_fine.ic(:,cn(1))>ic_threshold),...
                angl(hairc_fine.ic(:,cn(1))>ic_threshold,cn(1)),'r.');
            title(['Azimuth over time in ch. ' num2str(cn(1)) ' at cf = ' ...
                num2str(round(fc(cn(1)))) ' Hz'],'Fontsize',fontsize)
        case 4
            bar(-90:2:90,mean(h_all,2))
            hold on
            bar(-90:2:90,mean(h_ic,2),'r')
            title('Mean histogram of all fine-structure channels','Fontsize',fontsize)
            ymax = max(mean(h_all,2));
    end
    set(gca,'Fontsize',fontsize)
    if panel ~= 3
        set(gca,'XTick',s_pos)
        xlim([-93 93])
        ylim([0 ymax*1.1])
        xlabel('Azimuth [deg]','Fontsize',fontsize)
        ylabel('Frequency of occurence','Fontsize',fontsize)
        text (-80,ymax*.95,panellabel(panel),'Fontsize',fontsize+1,'FontWeight','bold')
    else
        set(gca,'YTick',s_pos)
        xlim([0 5.5])
        ylim([-95 95])
        ylabel('Azimuth [deg]','Fontsize',fontsize)
        xlabel('Time [s]','Fontsize',fontsize)
        rectangle('Position',[0.32,60,0.3,21],'FaceColor','white')
        text (0.38,69,panellabel(panel),'Fontsize',fontsize+1,'FontWeight','bold')
    end
end

example_output.angle_fine = angl;
example_output.IVS_fine = hairc_fine.ic;
example_output.histogram_angle_label = -90:2:90;
example_output.histograms_with_IVS = h_ic;
example_output.histograms_without_IVS = h_all;
