function definput=arg_baumgartner2014_pmv2ppp(definput)

definput.keyvals.p=ones(72,44);
definput.keyvals.rang=-90:5:269;
definput.keyvals.tang=[-30:5:70,80,100,110:5:210];
definput.keyvals.exptang=[];

definput.flags.print = {'noprint','print'};
definput.flags.chance = {'','chance'};
definput.flags.ppp = {'','QE_PE_EB','QE','PE','EB','absPE'};

end