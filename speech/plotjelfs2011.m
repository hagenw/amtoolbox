function op = plotjelfs2011(target_azim,database)
%PLOTJELFS2011  Plot output from jelf2011
%   Usage:  plotjelfs2011(target_azim,database);
%
%   `plotjelfs2011(target_azim,database)` will plot the output of the
%   `jelfs2011` binaural speech intelligibility advantage model for a
%   selected target azimuth angle. The masker position will move over a
%   full circle in the horizontal plane, and the output is visualized on
%   a polar plot.
%
%   The *database*  parameter selects the HRIR database to use. Please
%   see the help of |read_hrir|.
%
%   Examples:
%   ---------
%
%   Target angle of 0, database is `kemar`:::
%
%     plotjelfs2011(0,'kemar');
%
%   See also: jelfs2011, culling2005bmld, read_hrir
  
  step = 5;
  n_op = 360/step+1;
  op = zeros(n_op,3);
  angles = (0:step:360)'*pi/180;
  for direction_num = 1:n_op
    op(direction_num,:) = jelfs2011({target_azim, database}, ...
                                    {(direction_num-1)*step, database});
  end
  polar([angles angles angles], op);
  
end

function op = test_jelfs2011(targ_azim,int_azims,database)

  padding = zeros(1024,2);
  if strncmp(database,'siemens',7)
    fs = 48000;
  else
    fs = 44100;
  end
  
  int_ir = [];  
  target_ir = read_hrir(0,targ_azim,database); 
  
  for i=1:length(int_azims)
    int_ir = [int_ir; read_hrir(0,int_azims(i),database); padding];
  end
  
  int_ir = int_ir/sqrt(length(int_azims));
  
  op = jelfs2011([target_ir; padding],int_ir,fs);

end
