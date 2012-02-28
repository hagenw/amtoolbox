function op = do_prediction(targ_azim,int_azims,mode,database)
padding = zeros(1024,2);
if strncmp(database,'siemens',7)
    sf = 48000;
else
    sf = 44100;
end
int_ir = [];
target_ir = read_hrir(0,targ_azim,database); 
for i=1:length(int_azims)
  int_ir = [int_ir; read_hrir(0,int_azims(i),database); padding];
end
int_ir = int_ir/sqrt(length(int_azims));
op = jelfs2011([target_ir; padding],int_ir,sf,mode);
