function test_amtcache(varargin)

definput.import={'amtcache'}; % get the flags of amtcache
[flags,keyvals]  = ltfatarghelper({},definput,varargin);  % parse the input
% amtcache flags are in flags.cachemode now!

%% Cache a scalar
s=amtcache('get','scalar',flags.cachemode);
if isempty(s)
  disp('redo scalar')
  s=2;
  amtcache('set','scalar',s);
end
%% Cache a matrix
m=amtcache('get','matrix',flags.cachemode);
if isempty(m)
  disp('redo matrix')
  m=rand(10,10);
  amtcache('set','matrix',m);
end
%% Cache a structure
st=amtcache('get','struct',flags.cachemode);
if isempty(st)
  disp('redo struct')
  st.x=3;
  st.m=rand(10,100);
  st.string='sakdhfksjfhkjdsafk';
  amtcache('set','struct',st);
end

%% Cache multiple variables
[x,y,z]=amtcache('get','multi variables',flags.cachemode);
if isempty(x)
  disp('redo struct')
  x=3;
  y=rand(10,100);  
  z='sakdhfksjfhkjdsafk';
  amtcache('set','multi variables',x,y,z);
end

% amtcache('delete','fig6');
% amtcache('delete-full','amtcachetest/fig6');