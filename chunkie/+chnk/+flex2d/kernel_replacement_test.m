clear
clc


zk = 6; 
nu = 1/3; 

% t1 = 0;
% t2 = pi/10000;

coefs = [nu;0];

srcinfo = [];

srcinfo.r = [3; 0];
srcinfo.n = [1; 0];
srcinfo.d = [0; 1];



targinfo = [];


targinfo.r = [3; 1/1000];
targinfo.n = [1; 0];
targinfo.d = [0; 1];

fkern = @(s,t) chnk.helm2d.kern(zk, s, t, 'test kernel', coefs);
ans1 =  fkern(srcinfo, targinfo);

fkern1 = @(s,t) chnk.helm2d.kern(zk, s, t, 'free plate K21 first part', coefs);
ans2 = fkern1(srcinfo, targinfo);


abs(ans1-ans2)









