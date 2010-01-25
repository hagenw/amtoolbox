function outsig = lindemann(insig,show_band,inhibition_factor,weighting_factor,mono_processing,using_gaik)
% LINDEMANN Calculates a binaural activation pattern
%   Usage: outsig = lindemann()
%__________________________________________________________________________
%
% output = bimo(stimulus_file_name,show_band,inhibition_factor,weighting_factor,mono_processing,using_gaik)
%
% bimo calculates a binaural model for any given sound signal and plots the
% result.
%
% Options:
%   stimulus_file_name  - name of stimuli file in def.dir_stimuli to 
%                         process (without .mat)
%   show_band           - plot the results for this erb band
%   inhibition_factor   - constant inhibition factor (0.0 = no inhibition)
%   
%   weighting_factor    - Strength of the monaural processing of the
%                         Lindemann algorithm (see Lindemann 1986a, p. 1611
%                         eq. 9). It is weighting_factor < 1 and with
%                         weighting_factor = 0 no monaural processing is
%                         done at all.
%   mono_processing     - 1 or 0. Use also the mono signal for modeling.
%                         Note: this option only applies to the MEX-file
%                         implementation of the cross-correlation. In the
%                         m-file implementation of the Lindemann-model the
%                         monaural processing is automaticly applied using
%                         the monaural sensitivities with the
%                         weighting_factor. 
%   using_gaik          - 1 or 0. Use the Gaik (1993) algorithm to ad-
%                         vantage natural combinations of ITD and ILD 
%                         (measured with HRTFs)
%
% Output:
%   output              - A matrix containing the cross-correlation signal
%                         for every band and every time step n. The format
%                         of this matrix is output(n,m,band), where m
%                         denotes the correlation (delay line) time step.
%
% Example: bimo('Tonhalle',5,0.5,0.035,1,0);
%
% The steps of the binaural model to calculate its results are the
% following:
%
% 0) (optional)
% In a first step (if using_gaik == 1) naturally occuring binaural 
% parameters for a given location are calculated. This is done by 
% calculating the ILD and ITD for a given HRTF catalog. Natural binaural 
% parameters are of special interest, because they are handled different by
% the auditory system then unaturally occuring combinations of ITD and ILD 
% (see Gaik, 1993). 
% The results of these calculations are two level factors for the right and
% the left auditory channel. With these level factors the channels are
% multiplied calculating the cross-correlation between the channels. This
% multiplication is done in a way, that a naturally combination of an ITD
% and an ILD will get an ILD of zero in the Lindemann cross-correlation
% (and so a fusion of the lateralization image occure).
%
% 1)
% The given stimulus is filtered using a erb bank to get
% "def.erb_bank_end_band-def.erb_bank_start_band+1" frequency bands 
% containing a stimulus waveform.
%
% 2)
% For frequency bands with a high (> def.ild_start_band) center frequency
% the envelope is extracted and used for the calculation of the
% cross-correlation. This is done to simulate the decrease of phase locking
% of the auditory system for high frequencies (that means the temporal fine
% structur can only be used for low frequencies).
%
% 3)
% Half-wave rectification to simulate the correspondending behavior of the
% auditory hair cells.
%
% 4) 
% Calculation of the cross-correlation between the left and right channel.
% This is done using the model described in Lindemann (1986a) and Gaik
% (1993). These are extensions to the delay line model of Jeffres (1948).
% Lindemann has extended this model by a contralateral inhibition, which
% introduce the ILD to the model and explains the law of the first wave.
% Also monaural detectors were extended, to handle monaural signals (and
% some stimuli with a split off of the lateralization image).
% Gaik extended the Lindemann model by weighting of natural combinations
% of ITDs and ILDs. With this extension he could explain more experiments
% with ITDs and ILDs. Also he makes the inhibition independet of the signal
% amplitude.
% A detailed description of these cross-correlation steps is given in the
% binaural_correlation function.
% 
%
% ___Literature:___
%
% Gaik (1993):
%   "Combined evaluation of interaural time and intensity differences:
%   Psychoacoustic results and computer modeling", Werner Gaik, 1993 JASA
%
% Jeffres (1948):
%   "A place theory of sound localization", L.A. Jeffress, 1948 J. Comp. 
%   Physiol. Psychol.
%
% Lindemann (1986a):
%   "Extension of a binaural cross-correlation model by contralateral
%   inhibition. I. Simulation of lateralization for stationary signals", W.
%   Lindemann, 1986 JASA
%
% Lindemann (1986b):
%   "Extension of a binaural cross-correlation model by contralateral
%   inhibition. II. The law of the first wave front", W. Lindemann, 1986
%   JASA
%
% see also: bimo_config
%__________________________________________________________________________
%bi.mo                                                                v.0.3

%
% The binaural model (bimo) was first created by Wolfgang Hess at the
% Institue of Communication Acoustics using their long experience with
% binaural auditory models (see e.g. the first chapter of the book: Blauert
% "Communication Acoustics").
%
% Wolfgang Hess, 2005-jan
% Institute of Communication Acoustics, Ruhr-University Bochum, Germany
% Harman/Becker Automotive Systems, Karlsbad-Ittersbach, Germany
% wghess@gmx.de
%
% Hagen Wierstorf
% T-Labs, Berlin
% hagen.wierstorf@telekom.de
% 2009/05/28
%

%% NOTES
%
%
%

%% ------ Checking of input  parameters ---------------------------------------

error(nargchk(7,7,nargin));

if ~isnumeric(insig)
    error('%s: insig has to be numeric!',upper(mfilename));
end


%% ------ Variables -----------------------------------------------------------
fs = 44100;  % fs should be a parameter to the function
% Highest and lowest frequency to use for the erbfilterbank (this gives us 36 
% channels, channel 5-40)
flow = erbtofreq(5);
fhigh = erbtofreq(40); 

%% Make global variables/constants avaiable
% Read configuration struct (see bimo_config.m)
global def

% Store current directory and change to data directory
def.dir_current = pwd;
if ~strcmp(def.dir_current,def.dir_data)
    cd (def.dir_data)
end


%% Config section
len_y = 400; % 0,100,1000       % length of y axis. Why are it not calculatged from stimulus length?

% Start and stop of the elevation angle
% Is this ever used in this code?
%ele_start = 0;      
%ele_step = 5;       

mulfact = 10^4;     % Anhebung um 40dB, je größer mulfakt, desto zeitlich ausgedehnteres Muster 


% ------ Erb Bank -------------------------------------------------------------
% Generate an erb filterbank for simulation of the frequncy -> place
% transformation of the cochlea. This generates erb filterbank coefficients
% with a range from flow to fhigh.
% NOTE: Lindemann uses a bandpass filterbank after Duifhuis (1972) and
% Blauert and Cobben (1978). So for a copy of his model, this have to
% implemented as an option for the filterbank.
[b,a] = gammatone(erbspacebw(flow,fhigh),fs);
% Applying the erb filterbank to the signal
outsig = filterbank(b,a,insig);


%% Calculation of the level factors for the Gaik-algorithm
% This calculates natural binaural parameters (e.g. ITD and ILD) from a 
%>HRTF catalog in order to use the model described in Gaik (1993). In this
%>model combinations of ITDs (for envelope and fine structure) and ILDs 
%>were advantages, that are observable in a natural situation (using
%>HRTFs).
% In this section level factors are calculated using the measured ILDs and
%>ITDs for every frequency band. These level factors are later on applied 
%>at the calculation of the cross-correlation.
%if using_gaik == 1
if 0
    
    % step 1: Analysis of the binaural parameters of the HRTFs for every 
    %>direction

    % Progress bar
    progress_message = '=> Calculating natural combinations of ITD and ILD for the HRTFs';
    if def.show_progress_bar == 1
        h_waitbar = waitbar(0,progress_message);
    else
        fprintf(1,'%s:   000 / 100',progress_message);
    end

    % Calculate the correlation values of the HRTFs for single
    %>directions
    % Memory preallocation
    number_of_azimuths = length(def.azimuth_start:def.azimuth_step:(360-def.azimuth_step));
    hrtf_correlations.itd_indices = zeros(number_of_azimuths,def.erb_bank_end_band);
    hrtf_correlations.envelope_indices = zeros(number_of_azimuths,def.erb_bank_end_band);
    hrtf_correlations.ild_indices = zeros(number_of_azimuths,def.erb_bank_end_band);
    i = 1;
    for azimuth = def.azimuth_start:def.azimuth_step:(360-def.azimuth_step)

        % Progress bar
        if def.show_progress_bar == 1
            waitbar(azimuth/360,h_waitbar);
        else
            fprintf(1,'\b\b\b\b\b\b\b\b\b%3.0f / 100',azimuth/360*100);
        end

        % Spatialize the HRTF catalog. Choose only one direction (given
        %>by azimuth and elevation).
        spatialized_hrtf = spatialize(def.hrtf_catalog,[def.elevation;azimuth]);

        % Filter the spatialized HRTF signal with the erb filter bank
        %>and return an array with signals for every frequency band.
        filtered_hrtf = filter_bank(erb_bank,spatialized_hrtf);

        % Calculate binaural correlations (ILD, ITD, envelope
        %>correlation) for the HRTF signal
        [ itd_indices envelope_indices ild_indices ] = hrtf_analysis(filtered_hrtf);   

        % Store the calculated HRTF correlations
        hrtf_correlations.itd_indices(i,:) = itd_indices;
        hrtf_correlations.envelope_indices(i,:) = envelope_indices;
        hrtf_correlations.ild_indices(i,:) = ild_indices;

        % Increase counter
        i = i+1;
    end % of for azimuth
    clear i
    
    % step 2: estimation of the compensation parameters (level factors)
	gaik = prepare_gaik(hrtf_correlations);  
    
    % Progress bar
    if def.show_progress_bar == 1
        close(h_waitbar)
    else
        fprintf('\b\b\b\b\b\b\b\b\b%3.0f / 100\n',100);
    end
end
if using_gaik ~= 1
    % If we don't want to use the Gaik algorithm to advantage the natural
    %>combination of ITD and ILD, we have to set the level correction 
    %>factor to one (see Gaik 1993, p. 105/106).
    gaik.level_left = ones(def.number_of_itd_tabs,def.erb_bank_end_band);     
    gaik.level_right = ones(def.number_of_itd_tabs,def.erb_bank_end_band);     
end % of preparation


%% Calculate binaural activity patterns

% Progress bar
progress_message = '=> Preprocessing of the stimulus and calculating the cross-correlation';
if def.show_progress_bar == 1
    h_waitbar = waitbar(0,progress_message);
else
    fprintf(1,'%s:   000 / 100',progress_message);
end


% Extract the envelope, apply a half-wave rectification and calculate a
% running cross-correlation for every given frequency band
for band = def.erb_bank_start_band:def.erb_bank_end_band
	
    % ___E N V E L O P E   E X T R A C T I O N___
    % Only calculate the envelope for signals above the center_frequency
    %>given by the band "def.ild_start_band"
    if band > def.ild_start_band
        % Calculate the envelope using the hilbert transformation (Option
        %>'m2')
		signal_envelope = ([calculate_envelope(filtered_stimulus.signal_left(:,band),def.fs,100,'m2') calculate_envelope(filtered_stimulus.signal_right(:,band),def.fs,100,'m2')]);
	else 
		signal_envelope = ([filtered_stimulus.signal_left(:,band) filtered_stimulus.signal_right(:,band)]);
    end
    
    % NOTE: if the scaling of the data is enabled here, the spectral shape
    %>over the frequency bands are lost and you have to know which bands
    %>are of importance to you!!!
    % Scale signal_envelope (=> max(signal_envelope) = 0.999)
	%maximum = max(max(signal_envelope));
	%signal_envelope = signal_envelope*0.999./maximum;
    
    % ___H A L F - W A V E   R E C T I F I C A T I O N___
    % Half-wave rectification for both channels
    %>I don't know why we use eps and not 0.
    index = find(signal_envelope(:,1)<0); 
    signal_envelope(index,1) = eps;  
    index = find(signal_envelope(:,2)<0);
    signal_envelope(index,2) = eps;

    signal_left(:,band) = signal_envelope(:,1); 
    signal_right(:,band) = signal_envelope(:,2); 
    
    % ___C R O S S - C O R R E L A T I O N___ (m-file)
    % Calculate the cross-correlation after Lindemann (1986a) and Gaik
    %>(1993) using the m-file implementation (not the faster mex-file).
    if ~def.mex
        output(:,:,band) = binaural_correlation(signal_envelope,weighting_factor,inhibition_factor,gaik.level_left(:,band)',gaik.level_right(:,band)');
    end
    
    % Progress bar
    if def.show_progress_bar == 1
        waitbar((band-def.erb_bank_start_band)/(def.erb_bank_end_band-def.erb_bank_start_band),h_waitbar);
    else
        fprintf(1,'\b\b\b\b\b\b\b\b\b%3.0f / 100',(band-def.erb_bank_start_band)/(def.erb_bank_end_band-def.erb_bank_start_band)*100);
    end
    
end % of for band

% ___C R O S S _ C O R R E L A T I O N___ (mex-file)
% Calculate the cross-correlation after Lindemann (1986a) and Gaik (1993) 
%>using the mex-file implementation of the algorithm.
if def.mex
    output = mex_bimocorr(signal_right,signal_left,inhibition_factor,weighting_factor,gaik.level_right',gaik.level_left',mono_processing);
    % Resorting of output (really necessary?)
    output(1,1) = 0; %workaround
    temp = output;
    templen = floor(1000*length(signal_envelope)/def.fs);
    clear output
    % This can create an out of memory error for long (>1s) signals
    for i=1:def.erb_bank_end_band
        output(:,:,i) = temp(:,i:def.erb_bank_end_band:templen*def.erb_bank_end_band)';
    end
end

% Progress bar
if def.show_progress_bar == 1
    close(h_waitbar);
else
    fprintf('\b\b\b\b\b\b\b\b\b%3.0f / 100\n',100);
end


%% Plot XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

% extra file ? YES
% change structure, etc.

% This whole block calculates the length of the y axis
% (this influences primarly the time plot, but if we use 1000, then
% also the 3d plot is affected!)
if len_y == 0
    yend = 0;
    for i = 30:length(output(:,1))
        if sum(output(i,:),2) == 0
            yend = 1 + yend;
            if yend > 20
                ylength = i;
                break;
            end
        else    
            ylength = length(output(:,1));
            yend = 0;
        end
    end
elseif len_y == 1000
    sizeX = size(output,1);
    if sizeX < 1000
        output(sizeX+1:1000,1:181) = zeros(1000-sizeX,181); 
    end
    ylength = 1000;
elseif ((len_y <= 1000) && (len_y > 0))
    ylength = len_y;  
end


%% Animations ??    
if (def.generate_figure || def.generate_avi)
    
    % Name for the avi file
    outfile = ['G_' stimulus_file_name '_' int2str_(def.erb_bank_start_band,2,'0') 'to' int2str_(def.erb_bank_end_band,2,'0')...
               '_' int2str_(inhibition_factor*1000,3,'0') '_' int2str_(weighting_factor*1000,4,'0') '_' int2str_(mono_processing,1,'0') '_' int2str_(using_gaik,1,'0') '_000.avi'];
    % Check if the file already exists
    if exist(outfile,'file')
        % The following have to be completed so it will be able to handle a
        % random number of existing files.
        %outfile = regexprep(outfile,'000.avi','001.avi');_
        % Therefore we are deleting the file at the moment until file
        % handling works correctly.
        delete(outfile)
    end
    % Generate avi
    mov = avifile(outfile,'compression','none','quality',100,'fps',3);
    %
    map = gaik.itd;
    if def.average_num
        def.erb_bank_end_band = def.erb_bank_end_band-def.average_num;
    end
    for band=show_band:def.erb_bank_end_band
        
        if def.average_num
            initp = 1;
            for cband=band-def.average_num:band+def.average_num
                % disp(sprintf('Band %d \t',band));
                data = output(:,:,cband);
                for n=1:length(data(:,1))
                    rk(n,:)=data(n,:).*amp2db(mulfact.*max(data(n,:))+eps)./(max(data(n,:))+eps); % wg amp2db
                    index=find(rk(n,:)<0);                % Warum hier Halbwellengleichrichtung ???
                    rk(n,index)=0;
                    ycor(n,:)=rk(n,round(map(band,:)));
                end 
                rk=ycor+eps;
                shifto=-30;
                if initp == 1
                    returndata_n=rk./max(max(rk));
                    initp = 0;
                else
                    returndata_n=(returndata_n+rk./max(max(rk))).*0.5;
                end    
                clear rk;
            end
        else
            data = output(:,:,band);
            for n=1:length(data(:,1))
                rk(n,:)=data(n,:).*amp2db(mulfact.*max(data(n,:))+eps)./(max(data(n,:))+eps); % wg amp2db
                index=find(rk(n,:)<0);                
                rk(n,index)=0;
                ycor(n,:)=rk(n,round(map(band,:)));
            end 
            rk=ycor;
            shifto=-30;
            returndata_n=rk./max(max(rk));
            clear rk;
        end        
        to = length(returndata_n(:,1));
        %returndata_n = Interpolate2D(returndata_n,interpolate);
		figure('Position',[200 100 1224 868]);
		subplot('Position',[0.25 0.3 0.68 0.685])
	    %waterfall(returndata_n); colormap(gray)
        %surf(returndata_n);shading interp;colormap (gray);lighting phong;
		surf(returndata_n);shading interp;colormap (jet);lighting phong;
        grid on
        t = length(returndata_n(:,1));
        l = length(returndata_n(1,:));
   		set(gca,'XTick', [1 floor(1*l/12) floor(2*l/12) floor(3*l/12) floor(4*l/12) floor(5*l/12) floor(6*l/12) ...
                            floor(7*l/12) floor(8*l/12) floor(9*l/12) floor(10*l/12) floor(11*l/12) floor(l)]);
 		set(gca,'XtickLabel', ['-90';'-75';'-60';'-45';'-30';'-15';' 0 ';'15 ';'30 ';'45 ';'60 ';'75 ';'90 ']);
 		set(gca,'YTick', [1 floor(1*t/5)  floor(2*t/5) floor(3*t/5) floor(4*t/5) floor(5*t/5)]);
 		set(gca,'YtickLabel', 0:round(to/5):to);
		xlabel('lateralization[deg]','FontSize',14,'FontWeight','bold');
	    %ylabel('time [ms] =>','FontSize',10,'FontWeight','bold','rotation',82,'VerticalAlignment','bottom');
		ylabel('time[ms] =>','FontSize',12,'FontWeight','bold','VerticalAlignment','bottom','rotation',75);
		zlabel('activity =>','FontSize',12,'FontWeight','bold','VerticalAlignment','bottom');
		set(gca,'FontName','Times');
		axis([1 length(returndata_n(1,:)) 1 length(returndata_n(:,1))]);
        view(10,70)

		%Reverberation
		subplot('Position',[0.06 0.32 0.1 0.44])
		k_mean = max(returndata_n,[],2);
		k_mean_n = -k_mean;
		k_rf = 1:length(k_mean);
		cla;line('XData',k_mean_n,'YData',k_rf,'Color','r','LineWidth',2)
        %axis([min(k_mean_n) max(k_mean_n) 1 ylength])
        max_k_mean = 1.1;
  		axis([-max_k_mean 0 1 ylength])
 		set(gca,'XTick', [-max_k_mean -0.5*max_k_mean 0 ]);
 		set(gca,'XtickLabel', [max_k_mean; 0.5*max_k_mean; 0]);
 		set(gca,'YTick', [1 floor(1*t/5)  floor(2*t/5) floor(3*t/5) floor(4*t/5) floor(5*t/5)]);
 		set(gca,'YtickLabel', 0:round(to/5):to);
		xlabel('<=activity','FontSize',10,'FontWeight','bold','VerticalAlignment','top');
		ylabel('time[ms] =>','FontSize',10,'FontWeight','bold','VerticalAlignment','bottom');
		set(gca,'YAxisLocation','left','box','on');
        grid on;

		% activity
		subplot('Position',[0.25 0.05 0.62 0.1])
		k_mean = max(returndata_n,[],1);
		k_rf = 1:length(k_mean);
		cla;line('XData',k_rf,'YData',k_mean,'Color','b','LineWidth',2)
		%axis([1 length(k_mean(1,:)) min(k_mean) max(k_mean)])
        max_k_mean = 1; % 0.15
		axis([1 length(k_mean(1,:)) 0 1.1*max_k_mean])
   		set(gca,'XTick', [1 floor(1*l/12) floor(2*l/12) floor(3*l/12) floor(4*l/12) floor(5*l/12) floor(6*l/12) ...
                            floor(7*l/12) floor(8*l/12) floor(9*l/12) floor(10*l/12) floor(11*l/12) floor(l)]);
 		set(gca,'XtickLabel', ['-90';'-75';'-60';'-45';'-30';'-15';' 0 ';'15 ';'30 ';'45 ';'60 ';'75 ';'90 ']);
        %xlabel('lateralization[deg]','FontSize',10,'FontWeight','bold');
		ylabel('activity =>','FontSize',10,'FontWeight','bold');
		set(gca,'YAxisLocation','right','box','on');
        grid on

		text(-34+shifto,max_k_mean*1.5,stimulus_file_name,'color','b','FontWeight','bold','FontSize',10);
        % The following loads erb_bank.center_frequency, erb_bank.left_range and RightRange
		%temptext = [def.erb_bank_file,'.freq']; load(temptext,'-mat');
        if def.average_num
            temptext = ['band ', num2str(band-def.average_num), ' - ', num2str(band+def.average_num)];
		    text(-34+shifto,max_k_mean*1.0,temptext,'color','r','FontSize',10);
            
    		temptext = ['cntr frq ', num2str(round(erb_bank.center_frequency(band-def.average_num))),' - ', num2str(round(erb_bank.center_frequency(band+def.average_num))),' Hz'];
    		text(-34+shifto,max_k_mean*0.75,temptext,'color','r','FontSize',10);
        else
            temptext = ['band ', num2str(band)];
    		text(-34+shifto,max_k_mean*1.0,temptext,'color','r','FontSize',10);
    		temptext = ['cntr frq ', num2str(round(erb_bank.center_frequency(band))),' Hz'];
    		text(-34+shifto,max_k_mean*0.75,temptext,'color','r','FontSize',10);
        end        
		temptext = ['cs=',num2str(inhibition_factor),'  wf=',num2str(weighting_factor)];
		text(-34+shifto,max_k_mean*(0.5),temptext,'color','b','FontSize',10);
		if mono_processing == 0, temptext = ['monaural processors off']; else temptext = ['monaural processors on']; end
		text(-34+shifto,max_k_mean*(0.25),temptext,'color','b','FontSize',10);
		if using_gaik == 0, temptext = ['trading off']; else temptext = ['trading on']; end
		text(-34+shifto,max_k_mean*(0.0),temptext,'color','b','FontSize',10);
        
		if def.eps
            if def.average_num
		        eval(['print -depsc2 -tiff ' 'L_' stimulus_file_name '_' int2str_(def.azimuth_start,3,'0') '_' int2str_(def.elevation,3,'0') ...
					'_avg' int2str_(band-def.average_num,2,'0') '-' int2str_(band+def.average_num,2,'0') '_' int2str_(inhibition_factor*1000,3,'0') ...
					'_' int2str_(weighting_factor*1000,4,'0') '_' int2str_(mono_processing,1,'0') '_' int2str_(using_gaik,1,'0') '.eps']);
            else
  			    eval(['print -depsc2 -tiff ' 'L_' stimulus_file_name '_' int2str_(def.azimuth_start,3,'0') '_' int2str_(def.elevation,3,'0') ...
					'_' int2str_(band,2,'0') '_' int2str_(inhibition_factor*1000,3,'0') ...
					'_' int2str_(weighting_factor*1000,4,'0') '_' int2str_(mono_processing,1,'0') '_' int2str_(using_gaik,1,'0') '.eps']);
            end            
		end

        if def.png
            if def.average_num
		        eval(['print -dpng ' 'L_' stimulus_file_name '_' int2str_(def.azimuth_start,3,'0') '_' int2str_(def.elevation,3,'0') ...
					'_avg' int2str_(band-def.average_num,2,'0') '-' int2str_(band+def.average_num,2,'0') '_' int2str_(inhibition_factor*1000,3,'0') ...
					'_' int2str_(weighting_factor*1000,4,'0') '_' int2str_(mono_processing,1,'0') '_' int2str_(using_gaik,1,'0') '.png']);
            else
  			    eval(['print -dpng ' 'L_' stimulus_file_name '_' int2str_(def.azimuth_start,3,'0') '_' int2str_(def.elevation,3,'0') ...
					'_' int2str_(band,2,'0') '_' int2str_(inhibition_factor*1000,3,'0') ...
					'_' int2str_(weighting_factor*1000,4,'0') '_' int2str_(mono_processing,1,'0') '_' int2str_(using_gaik,1,'0') '.png']);
            end            
		end

        if def.generate_figure
            F = getframe(gcf);
   	        mov = addframe(mov,F);
            mov = close(mov);
            clear mov;
            [XI,YI] = meshgrid(-3:.125:3);            
            delete(outfile);
            break
        end
        set(gcf,'MenuBar','none');
        F = getframe(gcf);
        close(gcf)
   	    mov = addframe(mov,F);
    end %of for
    if exist('mov'), mov = close(mov); end;
end %of def.generate_figure | def.generate_avi

% IRs XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX 

if def.plot_signal 
    figure
    eval(['load ' stimulus_file_name]);

    for n=1:length(MAP(1,:))
        if (MAP(1,n) == 0)
            IRl = HRIR(:,n*2-1);
            IRr = HRIR(:,n*2);
            break
        end
    end % of for 

    rf=1:length(IRl);												
    rf=(rf-1)*(1/def.fs)*1000;  % 1000 == msec!!!
    [LMax,LMaxI] = max(abs(IRl));
    [RMax,RMaxI] = max(abs(IRr));
    LMaxI = LMaxI*1000/def.fs;
    RMaxI = RMaxI*1000/def.fs;
    LMaxIStr=['|max| @ ',num2str(LMaxI),' ms'];
    RMaxIStr=['|max| @ ',num2str(RMaxI),' ms'];
    
    SNRl = ceil(10*GetSNR(IRl,def.fs))/10;     
    SNRr = ceil(10*GetSNR(IRr,def.fs))/10;
    SNR = max(SNRl,SNRr);
    SNRlStr=['SNR = ',num2str(SNRl),' dB'];
    SNRrStr=['SNR = ',num2str(SNRr),' dB'];

    subplot('Position',[0.05 0.6 0.9 0.35])
    line('XData', rf(1:length(IRl)),'YData', IRl(1:length(IRl)),'Color', 'b');
    axis([0 ylength -1.1 1.1])
    xlabel('time [ms]','FontSize',12,'FontWeight','bold');
    text (1, 1.1,'Left channel','FontSize',14,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','bottom','Color','k')
    text (ylength/2, 1,LMaxIStr,'FontSize',12,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','top','Color','r')
    set(gca,'box','on');

    subplot('Position',[0.05 0.1 0.9 0.35])
    line('XData', rf(1:length(IRr)),'YData', IRr(1:length(IRr)),'Color', 'b');
    axis([0 ylength -1.1 1.1])
    xlabel('time [ms]','FontSize',12,'FontWeight','bold');
    text (1, 1.1,'Right channel','FontSize',14,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','bottom','Color','k')
    text (ylength/2, 1,RMaxIStr,'FontSize',12,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','top','Color','r')
    set(gca,'box','on');
    
    if def.eps
  	    eval(['print -depsc2 -tiff ' 'IR_'  stimulus_file_name '_' int2str_(def.azimuth_start,3,'0') '_' int2str_(def.elevation,3,'0') ...
                                 '_' int2str_(inhibition_factor*1000,3,'0') '_' int2str_(weighting_factor*1000,4,'0') '_' int2str_(mono_processing,1,'0') ...
                                 '_' int2str_(using_gaik,1,'0') '_' int2str_(ylength,4,'0') '.eps']);
    end
    
    if def.png
        eval(['print -dpng ' 'IR_' stimulus_file_name '_' int2str_(def.azimuth_start,3,'0') '_' int2str_(def.elevation,3,'0') ...
                                 '_' int2str_(inhibition_factor*1000,3,'0') '_' int2str_(weighting_factor*1000,4,'0') '_' int2str_(mono_processing,1,'0') ...
                                 '_' int2str_(using_gaik,1,'0') '_' int2str_(ylength,4,'0') '.png']);
    end
end

cd(def.dir_current);
return
