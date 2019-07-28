%TEST_ARCLENGTHSMOOTH
%
% demonstrate that the arclength function is not going to be smooth at the
% level of the discretization unless some rescaling of the panel is done 

iseed = 8675309;
rng(iseed);

addpaths_loc();

% define curve

cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 0;
pref = []; 
pref.k = 16;
narms = 5;
amp = 0.25;

% resolve curve adaptively

start = tic; chnkr = chunkfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);
fprintf('%5.2e s : time to build geo\n',t1)

% sort chunker 
chnkr = sort(chnkr);

dsdt = arclength(chnkr);
ipan = 1:3;
dsdt3 = dsdt(:,ipan); % arclength on 3 consecutive panels
x = lege.exps(chnkr.k)+1;

% plot arc length imagining panels to be length 2

tt = [x(:); 2+x(:); 4+x(:)];
figure(1)
subplot(1,2,1)
plot(tt,dsdt3(:))

% plot arc length scaling the panel lengths (h stores the known length of
% the panels in the underlying paramter space. to do this without that
% knowledge would probably require some resampling...)

h3 =chnkr.h(ipan);
h3sum = cumsum(h3);
tt = [x(:)*h3(1)/2; h3sum(1)+x(:)*h3(2)/2; h3sum(2)+x(:)*h3(3)/2];
subplot(1,2,2)
plot(tt,dsdt3(:))
