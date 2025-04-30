chunkerarcparamTest0();

function chunkerarcparamTest0()
% test arclength parameterization and reparamaterization


% get starting chnkr
cparams = []; 
cparams.maxchunklen = 0.3; 
cparams.eps = 1e-8;
fcurve = @(s) starfish(s,3,0.4);
chnkra = chunkerfunc(fcurve, cparams);

cparams = []; 
fcurve = @(s) starfish(s,4,0,[],pi,3);
chnkrb = chunkerfunc(fcurve, cparams);

chnkr = merge([chnkra,chnkrb]);

% get parameterization
param_data = chunkerarcparam_init(chnkr);
fcurvep = @(s) chunkerarcparam(s,param_data);

% check parameterization returns original points
ssa = arclengthfun(chnkra);
ssb = arclengthfun(chnkrb);
lena = sum(chnkra.wts(:));

ss = [ssa, ssb+lena];
rps= fcurvep(ss(:));

err_rs = norm(chnkr.r(:)  - rps(:));
fprintf('error in arclength r evaluation = %e\n\n',err_rs)
assert(err_rs< 1e-10)


% check that arclength parameterization is consistent
k = chnkr.k;
[xs,~,us,vs] = lege.exps(k);
dermat = (vs*[lege.derpol(us); zeros(1,k)]).';

a = 0.58; b = a + 0.2;

dermat=dermat*(2/(b-a));
ss = (b-a)*(xs+1)/2 +a;

[rp,dp,d2p] = fcurvep(ss);

fprintf('drds err = %e, drdds err = %e\n', ...
    norm(rp(:,:)*dermat - dp(:,:),1), norm(rp(:,:)*dermat*dermat - d2p(:,:),1))
fprintf('drdds err2 = %e, orthogonality err = %e\n\n',...
    norm(dp(:,:)*dermat - d2p(:,:),1), norm(sum(dp.*d2p,1),1))
assert(max([norm(rp(:,:)*dermat - dp(:,:),1),...
    norm(rp(:,:)*dermat*dermat - d2p(:,:),1),...
    norm(dp(:,:)*dermat - d2p(:,:),1), norm(sum(dp.*d2p,1),1)]) < 1e-8)




% check reparameterized chunker has correct area and length
[chnkr2,eps] = arc_param_chunker(chnkr);
area_err = abs(area(chnkr2) - area(chnkr));
len_err = abs(sum(chnkr2.wts(:))-sum(chnkr.wts(:)));
fprintf('reparameterize with eps= %e\n', eps)
fprintf('area error = %e, length error = %e\n\n',area_err,len_err)

assert(max(area_err,len_err) < 1e-8)

% plot coefficients
figure(2)
plot(1:16,log10(eps*ones(1,16)),'--','linewidth',2)
hold on
plot(log10(abs(param_data.cr(:,:))),'linewidth',1.5)
hold off
ylabel('log10 coefficient')
legend('parameterization tolerance')


% check that chnkr2 has an arclength parameterization

darc = arclengthdens(chnkr2);
lens = chunklen(chnkr2);
darc = darc./ (lens/2).';
fprintf('arclength param err = %e\n',norm(darc - 1,inf))
assert(norm(darc - 1,inf) < 1e-6)

end