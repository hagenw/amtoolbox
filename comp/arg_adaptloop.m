function definput=arg_adaptloop(definput)
 
  definput.keyvals.limit=10;
  definput.keyvals.minlvl=0;
  definput.keyvals.tau=[0.005 0.050 0.129 0.253 0.500];
  
  definput.groups.dau      = {'tau',[0.005 0.050 0.129 0.253 0.500]};
  definput.groups.breebaart = {'tau',linspace(0.005,0.5,5)};
  definput.groups.puschel  = {'tau',linspace(0.005,0.5,5),'limit',0};
