function definput=arg_drnl(definput)
 
  definput.keyvals.limit=10;
  definput.keyvals.minlvl=1e-5;
  definput.keyvals.tau=[0.005 0.050 0.129 0.253 0.500];
  
  definput.groups.dau     ={'limit',10,'minlvl',1e-5, ...
                      'tau',[0.005 0.050 0.129 0.253 0.500]};
  definput.groups.breebart={'limit',10,'minlvl',1e-5,...
                      'tau',linspace(0.005,0.5,5)};
  
