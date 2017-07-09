function definput=arg_lindemann1986_bincorr(definput)

definput.keyvals.c_s   = 0.3;
definput.keyvals.w_f   = 0.035;
definput.keyvals.M_f   = 6;
definput.keyvals.T_int = 5;
definput.keyvals.N_1   = 1;

definput.groups.stationary={'T_int',Inf,'N_1',17640'};
