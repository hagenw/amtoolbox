function data = data_baumgartner2014for2016(varargin)
%DATA_BAUMGARTNER2014FOR2016  Data from Baumgartner et al. (2014) for
%predictions using baumgartner2016.
%   Usage: data = data_baumgartner2014for2016(flag)
%
%   `data_baumgartner2014for2016(flag)` returns data from Baumgartner et al. (2014)
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
%     'redo'      DTFs and recalibrated sensitivity paramters of all listeners.
%
%   The fields in the output contains the following information
%
%     .id         listener ID
%
%     .S          listener-specific sensitivity parameter
%
%     .mrs        listener-specific task-induced response scatter (derived
%                 from central lateral response precision in baseline condition)
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
%     data_baumgartner2014for2016('pool');
%
%   To get all listener-specific data of the pool including experimental 
%   baseline data, use::
%
%     data_baumgartner2014for2016('baseline');
%
%   See also: baumgartner2014, exp_baumgartner2014

% AUTHOR : Robert Baumgartner

%% ------ Check input options --------------------------------------------

definput.keyvals.latseg = 0;
definput.keyvals.dlat = 30;

% Define input flags
% definput.flags.plot = {'noplot','plot'};
definput.flags.type = {'pool','baseline'};
% definput.flags.recalib = {'norecalib','recalib'};
definput.flags.HRTFformat = {'sofa','ari'};
definput.flags.recalib={'','recalib'};

definput.import={'baumgartner2016','amtcache'};

% Parse input options
[flags,kv]  = ltfatarghelper({'mrsmsp','gamma'},definput,varargin);

if flags.do_recalib
  flags.cachemode = 'redo';
end

% if flags.do_missingflag
%   flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
%              sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
%   error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
% end;
    

%% Listener pool (listener-specific SP-DTFs) 
if flags.do_pool || flags.do_baseline
  
  listeners = {'NH12';'NH15';'NH21';'NH22';'NH33';'NH39';'NH41';'NH42';... % NumChan
               'NH43';'NH46';'NH53';'NH55';'NH58';'NH62';'NH64';'NH68';...
               'NH71';'NH72';'NH14';'NH16';'NH17';'NH18';'NH57';...
               };
  data=cell2struct(listeners,'id',2);
             
    for ii = 1:length(data)
      
      data(ii).S = 0.5; % default sensitivity
      
      filename = fullfile(SOFAdbPath,'baumgartner2014',...
        ['ARI_' data(ii).id '_hrtf_M_dtf 256.sofa']);
      
      data(ii).Obj = SOFAload(filename);
      data(ii).fs = data(ii).Obj.Data.SamplingRate;
      
    end
  
    data = loadBaselineData(data,kv.latseg,kv.dlat);
  
    %% prior districution
    dang = 30; % angular width of segments
    Tmin = 5; % min. # of targets to estimate prior distribution
    edges = -90:dang:270;
    for ii = 1:length(data)
      try
        T = histcounts(data(ii).itemlist(:,6),edges);
        R = histcounts(data(ii).itemlist(:,8),edges);
      catch
        centers = edges(1:end-1)+diff(edges)/2;
        T = hist(data(ii).itemlist(:,6),centers);
        R = hist(data(ii).itemlist(:,8),centers);
      end
      T(T<Tmin) = nan;
      RvT = R./T;
      RvT(isnan(RvT)) = 1;
      data(ii).priordist.y = RvT;
      data(ii).priordist.x = edges(1:end-1)+dang/2;
    end
    
  %% Calibration of S
  
  % Define cache name according to settings for auditory periphery model

  cachename = ['calibration_g' num2str(kv.gamma,'%u') ...
      '_mrs' num2str(kv.mrsmsp,'%u') ...
      '_do' num2str(kv.do,'%u') ...
      '_tar' num2str(kv.SPL,'%u') 'dB_tem' num2str(kv.SPLtem,'%u') 'dB_'...
      flags.fbank];
  if flags.do_gammatone
    cachename = [cachename '_'  num2str(1/kv.space,'%u') 'bpERB'];
    if flags.do_middleear; cachename = [cachename '_middleear']; end
    if flags.do_ihc; cachename = [cachename '_ihc']; end
  else % zilany
    cachename = [cachename '_' flags.fibertypeseparation];
  end
  if kv.prior > 0 
    cachename = [cachename '_prior' num2str(kv.prior,'%u')];
  end
  if kv.tiwin < 0.5
    cachename = [cachename '_tiwin' num2str(kv.tiwin*1e3) 'ms']; 
  end
  cachename = [cachename '_mgs' num2str(kv.mgs)]; 
  
  c = amtcache('get',cachename,flags.cachemode);
  if isempty(c) || not(isequal(c.kv,kv))
%   if not(exist('baumgartner2014for2016calibration.mat','file')) || flags.do_recalib
    
%     data = loadBaselineData(data,kv.latseg,kv.dlat);
    
    % reset listener-specific MRS to default
    for ii = 1:length(data)
      data(ii).mrs = kv.mrsmsp;
    end
    
    amtdisp('Calibration procedure started. Please wait!','progress')
%     data = baumgartner2016calibration(data,kv,flags);
    data = baumgartner2016calibration(data,'latseg',kv.latseg,...
      'SPLtem',kv.SPLtem,'gamma',kv.gamma,'prior',kv.prior,...
      flags.fbank,flags.fibertypeseparation,'tiwin',kv.tiwin,'mgs',kv.mgs);
    
%     data_all = data;
%     data = rmfield(data,{'Obj','itemlist','fs','target','response'}); % reduce filesize
%     save(fullfile(amtbasepath,'modelstages','baumgartner2014for2016calibration.mat'),'data')
%     data = data_all;
    
    c.data = rmfield(data,{'Obj','fs','itemlist','target','response'}); % reduce filesize
    c.kv = kv;
    amtcache('set',cachename,c)
    
  else
    
%     if flags.do_baseline
%       data = loadBaselineData(data,kv.latseg,kv.dlat);
%     end
    
%     c = load('baumgartner2014for2016calibration.mat');
      
    for ss = 1:length(data)
      for ii = 1:length(c.data)
        if strcmp(data(ss).id,c.data(ii).id)
          data(ss).S = c.data(ii).S;
          data(ss).mrs = c.data(ii).mrs;
          if isfield(c.data,'prior')
            data(ss).prior = c.data(ii).prior;
          else
            data(ss).prior = kv.prior;
          end
        end
      end
    end
    
  end 

end
    


end



function s = loadBaselineData(s,latseg,dlat)

% latseg = 0;%[-20,0,20]; 
% dlat = 30;%10;

% Experimental baseline data
numchan = data_goupell2010('BB');
methods = data_majdak2010('Learn_M');
spatstrat = data_majdak2013('BB');
ctcL = data_majdak2013ctc('Learn');

for ll = 1:length(s)
  
  s(ll).itemlist = [];
  
  s(ll).itemlist = [s(ll).itemlist ; numchan(ismember({numchan.id},s(ll).id)).mtx];
  s(ll).itemlist = [s(ll).itemlist ; methods(ismember({methods.id},s(ll).id)).mtx];
  s(ll).itemlist = [s(ll).itemlist ; spatstrat(ismember({spatstrat.id},s(ll).id)).mtx];
  s(ll).itemlist = [s(ll).itemlist ; ctcL(ismember({ctcL.id},s(ll).id)).mtx];
  
  s(ll).pe_exp = localizationerror(s(ll).itemlist,'rmsPmedianlocal');
  s(ll).qe_exp = localizationerror(s(ll).itemlist,'querrMiddlebrooks');
  s(ll).mrs = localizationerror(s(ll).itemlist,'precLcentral');
  
  for ii = 1:length(latseg)
    
    latresp = s(ll).itemlist(:,7);
    idlat = latresp <= latseg(ii)+dlat & latresp > latseg(ii)-dlat;
    mm2 = s(ll).itemlist(idlat,:);
    
    s(ll).pe_exp_lat(ii) = localizationerror(mm2,'rmsPmedianlocal');
    s(ll).qe_exp_lat(ii) = localizationerror(mm2,'querrMiddlebrooks');
    
    s(ll).target{ii} = mm2(:,6); % polar angle of target
    s(ll).response{ii} = mm2(:,8); % polar angle of response
    s(ll).Ntar{ii} = length(s(ll).target{ii});

  end

  
end

end