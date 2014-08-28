function data = data_baumgartner2014(varargin)
%DATA_BAUMGARTNER2014  Data from Baumgartner et al. (2014)
%   Usage: data = data_baumgartner2014(flag)
%
%   `data_baumgartner2014(flag)` returns data from Baumgartner et al. (2014)
%   describing a model for sound localization in sagittal planes (SPs)
%   on the basis of listener-specific directional transfer functions (DTFs).
%
%   The flag may be one of:
%
%     'pool'      DTFs and calibration data of the pool. This is the
%                 default.
%
%     'baseline'  Same as 'pool', but also with experimental data for
%                 baseline condition.
%
%     'recalib'   DTFs and recalibrated sensitivity paramters of all listeners.
%
%   The fields in the output contains the following information
%
%     .id         listener ID
%
%     .S          listener-specific sensitivity parameter
%
%     .Obj        DTF data in SOFA Format
%
%     .pe_exp     experimental local polar RMS error
%
%     .qe_exp     experimental quadrant error rate
%
%     .target     experimental target angles
%
%     .response   experimental response angles
%
%     .itemlist   experimental item list. Columns denote:
%                 1:4 ... azi_target,ele_target,azi_response,ele_response
%                 5:8 ... lat_target,pol_target,lat_response,pol_response
%                 9   ... F/B-Confusion resolved pol_response
%
%   Requirements: 
%   -------------
%
%   1) SOFA API from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
% 
%   2) Data in hrtf/baumgartner2014
%
%   Examples:
%   ---------
%
%   To get all listener-specific data of the pool, use::
%
%     data_baumgartner2014('pool');
%
%   To get all listener-specific data of the pool including experimental 
%   baseline data, use::
%
%     data_baumgartner2014('baseline');
%
%   See also: baumgartner2014, exp_baumgartner2014

% AUTHOR : Robert Baumgartner

%% ------ Check input options --------------------------------------------

% Define input flags
% definput.flags.plot = {'noplot','plot'};
definput.flags.type = {'pool','baseline'};
definput.flags.recalib = {'norecalib','recalib'};
definput.flags.HRTFformat = {'sofa','ari'};

definput.import={'baumgartner2014'};

% Parse input options
[flags,kv]  = ltfatarghelper({'mrsmsp','gamma'},definput,varargin);

% if flags.do_missingflag
%   flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
%              sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
%   error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
% end;
    

%% Listener pool (listener-specific SP-DTFs) 
if flags.do_pool || flags.do_baseline
  
  listeners = {'NH12';'NH15';'NH21';'NH22';'NH33';'NH39';'NH41';'NH42';'NH43';...
               'NH46';'NH53';'NH55';'NH58';'NH62';'NH64';'NH68';'NH71';'NH72';...
               'NH14';'NH16';'NH17';'NH18';'NH57'};
  data=cell2struct(listeners,'id',2);
             
    for ii = 1:length(data)
      
      data(ii).S = 0.5; % default sensitivity
      
      filename = fullfile(SOFAdbPath,'baumgartner2014',...
        ['ARI_' data(ii).id '_hrtf_M_dtf 256.sofa']);
      
      data(ii).Obj = SOFAload(filename);
      data(ii).fs = data(ii).Obj.Data.SamplingRate;
      
    end
  
  
  %% Calibration of S
  if not(exist('baumgartner2014calibration.mat','file')) || flags.do_recalib
    
    data = loadBaselineData(data);
    disp('Calibration procedure started. Please wait!')
    data = baumgartner2014calibration(data,kv);
    
    data_all = data;
    data = rmfield(data,{'Obj','itemlist','fs','target','response'}); % reduce filesize
    save(fullfile(amtbasepath,'modelstages','baumgartner2014calibration.mat'),'data')
    data = data_all;
    
  else
    
    if flags.do_baseline
      data = loadBaselineData(data);
    end
    
    c = load('baumgartner2014calibration.mat');
      
    for ss = 1:length(data)
      for ii = 1:length(c.data)
        if strcmp(data(ss).id,c.data(ii).id)
          data(ss).S = c.data(ii).S;
        end
      end
    end
    
  end 

end
    


end



function s = loadBaselineData(s)

latseg = [-20,0,20]; 
dlat = 10;

% Experimental baseline data
numchan = data_goupell2010('BB');
methods = data_majdak2010('Learn_M');
spatstrat = data_majdak2013('BB');
% ctc = data_majdak2013ctc('A');
ctcL = data_majdak2013ctc('Learn');

for ll = 1:length(s)
  
  s(ll).itemlist = [];
  
  s(ll).itemlist = [s(ll).itemlist ; numchan(ismember({numchan.id},s(ll).id)).mtx];
  s(ll).itemlist = [s(ll).itemlist ; methods(ismember({methods.id},s(ll).id)).mtx];
  s(ll).itemlist = [s(ll).itemlist ; spatstrat(ismember({spatstrat.id},s(ll).id)).mtx];
%   s(ll).itemlist = [s(ll).itemlist ; ctcA(ismember({ctcA.id},s(ll).id)).mtx];
%   s(ll).itemlist = [s(ll).itemlist ; ctcB(ismember({ctcB.id},s(ll).id)).mtx];
  s(ll).itemlist = [s(ll).itemlist ; ctcL(ismember({ctcL.id},s(ll).id)).mtx];
  
  s(ll).pe_exp = localizationerror(s(ll).itemlist,'rmsPmedianlocal');
  s(ll).qe_exp = localizationerror(s(ll).itemlist,'querrMiddlebrooks');   
  
  for ii = 1:length(latseg)
    
    latresp = s(ll).itemlist(:,7);
    idlat = latresp <= latseg(ii)+dlat & latresp > latseg(ii)-dlat;
    mm2 = s(ll).itemlist(idlat,:);
    
    s(ll).pe_exp_lat(ii) = localizationerror(mm2,'rmsPmedianlocal');
    s(ll).qe_exp_lat(ii) = localizationerror(mm2,'querrMiddlebrooks');
    
    s(ll).target{ii} = mm2(:,6); % polar angle of target
    s(ll).response{ii} = mm2(:,8); % polar angle of response
    s(ll).Ntargets{ii} = length(s(ll).target{ii});

  end

  
end

end