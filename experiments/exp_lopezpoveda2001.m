function exp_lopezpoveda2001(varargin)
%LOPEZPOVEDA2001   Iso-itensity curves of the DRNL Lopez-Poveda and Meddis (2001)
%   Usage: lopezpoveda2001_fig3bc;
%
%   Isointensity response of the linear, nonlinear and summed response of the
%   DRNL filter
%
%   This script reproduces figure 3b and c of Lopez-Poveda and Meddis 2001
%
%R  lopezpoveda2001hnc
  
%  AUTHOR: Katharina Egger

%% ------ Check input options --------------------------------------------

  definput.flags.type = {'fig3b','fig3c'};
  definput.flags.plot = {'noplot','plot'};

  % Parse input options
  [flags,keyvals]  = ltfatarghelper({},definput,varargin);
        
%% parameter set of YO, table I at 1000Hz
  exppars= {...
      'flow',1000,...
      'fhigh',1000,...
      'basef',1000,...
      'bwmul',1,...      
      'lin_fc', [log10(945) 0],...
      'lin_bw', [log10(240) 0],...
      'lin_gain', [log10(520) 0],...
      'lin_lp_cutoff', [log10(945) 0],...      
      'nlin_ngt_before', 3,...
      'nlin_ngt_after', 3,...
      'nlin_nlp', 3,...
      'nlin_fc_before', [log10(1000) 0],...
      'nlin_fc_after', [log10(1000) 0],...
      'nlin_bw_before', [log10(175) 0],...
      'nlin_bw_after', [log10(175) 0],...
      'nlin_lp_cutoff', [log10(1000) 0],...      
      'nlin_a', [log10(4598) 0],...
      'nlin_b', [log10(0.130) 0],...
      'nlin_c', [log10(0.25) 0] };
  
%% input signal
% 50ms pure tones, sampled at 100kHz
fs = 100000;
T = 0.05;       
t = (0:1/fs:T - 1/fs)';
fsig = 250:25:1750;    

if flags.do_fig3b
  %% Lopez-Poveda and Meddis 2001, Figure 3, b)
  level = 10^(30/20);
  
  result = zeros(1,length(fsig));
  lin = zeros(1,length(fsig));
  nlin = zeros(1,length(fsig));
  
  for ii = 1:length(fsig)
    
    insig = sin(2*pi*fsig(ii).*t)*(2^0.5) * level;  
    
    hp_fir = headphonefilter(fs);
    insig = filter(hp_fir,1,insig);
    
    [y_lin, ~] = drnl(insig, fs, exppars{:},'linonly');
    
    [y_nlin, ~] = drnl(insig, fs, exppars{:},'nlinonly');
    
    y_lin = gaindb(y_lin,-50);   % undo the 50dB gain in the drnl (not used in Lopez-Poveda and Meddis 2001)
    y_nlin = gaindb(y_nlin,-50);
    outsig = y_lin + y_nlin;
    
    result(1,ii) = rms(outsig(floor(length(insig)/2):end));
    lin(1,ii) = rms(y_lin(floor(length(insig)/2):end));
    nlin(1,ii) = rms(y_nlin(floor(length(insig)/2):end));
  end
  
  plot(fsig,result)
  hold on
  plot(fsig,lin,'-.g')
  plot(fsig,nlin,':r')
  set(gca,'YScale','log')
  % grid on
  set(gca,'XLim',[250 1750],'Layer','top')
  set(gca,'YLim',[1e-07 1e-03],'Layer','top')
  title('30 dB SPL')
  xlabel('Frequency (Hz)')
  ylabel('DRNL filter output (m/s)')
  
end;


if flags.do_fig3c
  %% Lopez-Poveda and Meddis 2001, Figure 3, c)

  result = zeros(1,length(fsig));
  lin = zeros(1,length(fsig));
  nlin = zeros(1,length(fsig));
  
  level = 10^(85/20);
  
  for ii = 1:length(fsig)
    
    insig = sin(2*pi*fsig(ii).*t)*(2^0.5)* level;
    hp_fir = headphonefilter(fs);
    insig = filter(hp_fir,1,insig);
    
    [y_lin, ~] = drnl(insig, fs, exppars{:},'linonly');
    
    [y_nlin, ~] = drnl(insig, fs, exppars{:},'nlinonly');

    y_lin = gaindb(y_lin,-50);
    y_nlin = gaindb(y_nlin,-50);
    outsig = y_lin + y_nlin;
        
    result(1,ii) = rms(outsig(floor(length(insig)/2):end));
    lin(1,ii) = rms(y_lin(floor(length(insig)/2):end));
    nlin(1,ii) = rms(y_nlin(floor(length(insig)/2):end));
    
  end
  
  plot(fsig,result)
  hold on
  plot(fsig,lin,'-.g')
  plot(fsig,nlin,':r')
  set(gca,'YScale','log')
  % grid on
  set(gca,'XLim',[250 1750],'Layer','top')
  set(gca,'YLim',[1e-05 1e-01],'Layer','top')
  title('85 dB SPL')
  xlabel('Frequency (Hz)')
  ylabel('DRNL filter output (m/s)')
  legend('DRNL output', 'Linear path output', 'Nonlinear path output')
  %set(gcf,'Position',[50,50,270,760])
  %set(gcf,'Position',[50,50,400,760])
  
end;