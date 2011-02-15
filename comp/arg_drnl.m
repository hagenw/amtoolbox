function definput=arg_drnl(definput)
 
  definput.flags.middleear={'middleear','nomiddleear'};

  definput.keyvals.flow=80;
  definput.keyvals.fhigh=8000;
  definput.keyvals.basef=[];

  definput.keyvals.lin_ngt = 2; 
  definput.keyvals.lin_nlp = 4; 
  definput.keyvals.lin_fc = [-0.06762 1.01679];
  definput.keyvals.lin_bw = [  .03728 .78563];
  definput.keyvals.lin_gain = [4.20405 -.47909];
  definput.keyvals.lin_lp_cutoff = [.06762 1.01679 ];
  
  definput.keyvals.nlin_ngt_before = 2;
  definput.keyvals.nlin_ngt_after = 2;
  definput.keyvals.nlin_nlp = 1;
  definput.keyvals.nlin_fc_before = [-0.05252 1.01650];
  definput.keyvals.nlin_fc_after  = [-0.05252 1.01650 ];
  definput.keyvals.nlin_bw_before = [-0.03193 .77426 ];
  definput.keyvals.nlin_bw_after  = [-0.03193 .77426 ];
  definput.keyvals.nlin_lp_cutoff = [-0.05252 1.01 ];
  
  definput.keyvals.nlin_a = [1.40298 .81916 ];
  definput.keyvals.nlin_b = [1.61912 -.81867 ];
  definput.keyvals.nlin_c = [-.60206 0];
  definput.keyvals.nlin_d = 1;

  definput.keyvals.compresslimit = [];
  
  definput.groups.jepsen2008={...
      'lin_bw',          [ .03728   .75 ],...
      'lin_lp_cutoff',   [-0.06762 1.01 ],...  
      'nlin_bw_before' , [-0.03193  .7  ],...
      'nlin_bw_after'  , [-0.03193  .7  ],...
      'compresslimit', 1500
      };
