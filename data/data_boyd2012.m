function data = data_boyd2012
%DATA_BOYD2012 - Data from Boyd et al. (JASA-EL, 2012)
%   Usage: data = data_boyd2012
%
%   Mean externalization scores of NH listeners extracted from top panels 
%   (1 talker condition) of Fig. 1 
%
%   Output parameters:
%     data    : structure with fields
%                 ID ... subject ID
%                 Resp ... externalization responses 
%                 BRIR ... binaural room impulse responses for 4 positions
%
%   References: boyd2012

% AUTHOR: Robert Baumgartner, Acoustics Research Institute, Vienna, Austria

% Data extracted from Fig. 1 (using WebPlotDigitizer)
% data.reference = 100;
% data.mix = 100:-25:0;
% data.ITE.BB = [100,95,58,45,33];
% data.BTE.BB = [63,66,38,38,24];
% data.ITE.LP = [72,75,42,26,12];
% data.BTE.LP = [61,62,33,24,16];

fn = amt_load('boyd2012','data.xlsx');
num = xlsread(fn,1,'B5:I52');
for ii = 1:7
  data(ii).ID = ['00',num2str(ii)];
  data(ii).Resp.mix = [nan,100:-25:0];
  rows = (1:6)+(ii-1)*7;
  if exist(fullfile(amt_load('boyd2012',''),data(ii).ID,'imp1.wav'),'file')
    for ipos = 4:-1:1
      ifn = fullfile(amt_load('boyd2012',''),data(ii).ID,['imp',num2str(ipos)]);
      [data(ii).BRIR.ITE(:,:,ipos),data(ii).BRIR.fs] = audioread([ifn,'.wav']);
      data(ii).BRIR.BTE(:,:,ipos) = audioread([ifn,'BTE.wav']);
      data(ii).BRIR.noHead(:,:,ipos) = audioread([ifn,'NH.wav']);
    end
  end
  data(ii).Resp.ITE_BB_oneT = num(rows,1);
  data(ii).Resp.ITE_BB_oneT = num(rows,1);
  data(ii).Resp.ITE_BB_fourT = num(rows,2);
  data(ii).Resp.BTE_BB_oneT = num(rows,3);
  data(ii).Resp.BTE_BB_fourT = num(rows,4);
  data(ii).Resp.ITE_LP_oneT = num(rows,5);
  data(ii).Resp.ITE_LP_fourT = num(rows,6);
  data(ii).Resp.BTE_LP_oneT = num(rows,7);
  data(ii).Resp.BTE_LP_fourT = num(rows,8);
end

end