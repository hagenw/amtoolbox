function definput = arg_baumgartner2016calibration( definput )

definput.keyvals.TolX = 0.005;
definput.keyvals.MaxIter = 100;
definput.keyvals.Srange = [-1,2.5];%[eps,10];
definput.keyvals.prange = [0,1];
definput.keyvals.latseg = 0;
definput.keyvals.dlat = 30;
definput.keyvals.c = {};

definput.flags.recalib={'','recalib'};
definput.flags.prior = {'','calibprior'};
definput.flags.optimization = {'fminbnd','fminsearch','search'};

end

