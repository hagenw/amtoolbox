%DEMO_JELFS2011  Binaural speech intelligibility advantage
%
%   `demo_jelfs2011` will plot the output of the
%   `jelfs2011` binaural speech intelligibility advantage model for a
%   target azimuth angle of 0 deg. The masker position will move over a
%   full circle in the horizontal plane, and the output is visualized on
%   a polar plot. The KEMAR HRTF dataset is used.
%
%
%   See also: jelfs2011, culling2005bmld

target_azim=0;
database='kemar';

step = 5;
n_op = 360/step+1;
op = zeros(n_op,3);
angles = (0:step:360)'*pi/180;
for direction_num = 1:n_op
  op(direction_num,:) = jelfs2011({target_azim, database}, ...
                                  {(direction_num-1)*step, database});
end
polar([angles angles angles], op);