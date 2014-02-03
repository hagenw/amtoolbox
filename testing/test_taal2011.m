f=greasy;
f=[greasy;greasy];

L=numel(f);

x=f;
y=setdbspl(randn(L,1),dbspl(f)-20);

d=taal2011(x,y,16000);
d_ref=ref_stoi(x,y,16000);
d_ref_1=ref_stoi_1(x,y,16000);

% Exact reference
(d_ref-d_ref_1)/d_ref

% Modified with respect to the reference
(d-d_ref)/d_ref

% Test for multi-signal
x=randn(20000,2);
y=x+0.1*randn(20000,2);

d_ref(1)= ref_stoi_1(x(:,1),y(:,1),10000);
d_ref(2)= ref_stoi_1(x(:,2),y(:,2),10000);
d       = taal2011(x,y,10000);

(d-d_ref)./d_ref