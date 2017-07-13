function test_amt_cache(varargin)

definput.import={'amt_cache'}; % get the flags of amt_cache
[flags,keyvals]  = ltfatarghelper({},definput,varargin);  % parse the input
% amt_cache flags are in flags.cachemode now!

%% Cache a scalar
s=amt_cache('get','scalar',flags.cachemode);
if isempty(s)
  disp('redo scalar')
  s=2;
  amt_cache('set','scalar',s);
end
%% Cache a matrix
m=amt_cache('get','matrix',flags.cachemode);
if isempty(m)
  disp('redo matrix')
  m=rand(10,10);
  amt_cache('set','matrix',m);
end
%% Cache a structure
st=amt_cache('get','struct',flags.cachemode);
if isempty(st)
  disp('redo struct')
  st.x=3;
  st.m=rand(10,100);
  st.string='sakdhfksjfhkjdsafk';
  amt_cache('set','struct',st);
end

%% Cache multiple variables
[x,y,z]=amt_cache('get','multi variables',flags.cachemode);
if isempty(x)
  disp('redo struct')
  x=3;
  y=rand(10,100);  
  z='sakdhfksjfhkjdsafk';
  amt_cache('set','multi variables',x,y,z);
end

% amt_cache('delete','fig6');
% amt_cache('delete-full','amt_cachetest/fig6');