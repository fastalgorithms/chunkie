flagselfTest0();


function flagselfTest0()
% FLAGSELFTEST test the routine for flagging when sources and targets
% overlap


nsrc = 40000;
xtarg = linspace(0,1,floor(sqrt(nsrc))); 
ytarg = linspace(0,2,floor(sqrt(nsrc)));
[xxsrc,yysrc] = meshgrid(xtarg,ytarg);
srcs = zeros(2,length(xxsrc(:)));
srcs(1,:) = xxsrc(:); srcs(2,:) = yysrc(:);

nsrc = size(srcs,2);


targs = [srcs, [1;2].*rand(2,100)];
P = randperm(size(targs,2));
targs = targs(:,P)+1e-15*(2*rand(size(targs))-1);

tic
flagslf = chnk.flagself(srcs, targs);
toc
fprintf('number of flagged points: %d\n', size(flagslf,2))
assert(size(flagslf,2) == nsrc)


errs = sum(vecnorm(srcs(:,flagslf(1,:)) - targs(:,flagslf(2,:)))>1e-10);
fprintf('number of incorrectly flagged points: %d\n', errs)
assert(norm(srcs(:,flagslf(1,:)) - srcs(:,P(flagslf(2,:))),1)<1e-6)

end


