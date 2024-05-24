%test_helmdiffgreen test the difference kernel functions 
% using finite differences 
% 

clearvars; clc
rng(1234);

ns = 3;
nt = 4;
k = 1.3;

src0 = randn(2,ns);
targ0 = randn(2,nt);

% test derivatives

[val0,grad0,hess0,der30,der40] = chnk.helm2d.helmdiffgreen(k,src0,targ0);

for j = 1:5
    h = 10^(-j);
    dx = h*[1;0];
    targ1 = targ0 + dx;
    [val1,grad1,hess1,der31,der41] = chnk.helm2d.helmdiffgreen(k,src0,targ1);

    errdx = norm(ones(size(val0)) - 2*(val1-val0)/h./(grad0(:,:,1)+grad1(:,:,1)));
    errdxx = norm(ones(size(val0)) - 2*(grad1(:,:,1)-grad0(:,:,1))/h./(hess0(:,:,1)+hess1(:,:,1)));
    errdxy = norm(ones(size(val0)) - 2*(grad1(:,:,2)-grad0(:,:,2))/h./(hess0(:,:,2)+hess1(:,:,2)));
    errdxxx = norm(ones(size(val0)) - 2*(hess1(:,:,1)-hess0(:,:,1))/h./(der30(:,:,1)+der31(:,:,1)));
    errdxxy = norm(ones(size(val0)) - 2*(hess1(:,:,2)-hess0(:,:,2))/h./(der30(:,:,2)+der31(:,:,2)));
    errdxyy = norm(ones(size(val0)) - 2*(hess1(:,:,3)-hess0(:,:,3))/h./(der30(:,:,3)+der31(:,:,3)));
    errdxxxx = norm(ones(size(val0)) - 2*(der31(:,:,1)-der30(:,:,1))/h./(der40(:,:,1)+der41(:,:,1)));
    errdxxxy = norm(ones(size(val0)) - 2*(der31(:,:,2)-der30(:,:,2))/h./(der40(:,:,2)+der41(:,:,2)));
    errdxxyy = norm(ones(size(val0)) - 2*(der31(:,:,3)-der30(:,:,3))/h./(der40(:,:,3)+der41(:,:,3)));
    errdxyyy = norm(ones(size(val0)) - 2*(der31(:,:,4)-der30(:,:,4))/h./(der40(:,:,4)+der41(:,:,4)));
    
    dx = h*[0;1];
    targ1 = targ0 + dx;
    [val1,grad1,hess1,der31,der41] = chnk.helm2d.helmdiffgreen(k,src0,targ1);

    errdy = norm(ones(size(val0)) - 2*(val1-val0)/h./(grad0(:,:,2)+grad1(:,:,2)));
    errdyy = norm(ones(size(val0)) - 2*(grad1(:,:,2)-grad0(:,:,2))/h./(hess0(:,:,3)+hess1(:,:,3)));
    errdyyy = norm(ones(size(val0)) - 2*(hess1(:,:,3)-hess0(:,:,3))/h./(der30(:,:,4)+der31(:,:,4)));
    errdyyyy = norm(ones(size(val0)) - 2*(der31(:,:,4)-der30(:,:,4))/h./(der40(:,:,5)+der41(:,:,5)));

    fprintf('%5.2e : err in dx\n',errdx)    
    fprintf('%5.2e : err in dy\n',errdy)    
    fprintf('%5.2e : err in dxx\n',errdxx)    
    fprintf('%5.2e : err in dxy\n',errdxy)    
    fprintf('%5.2e : err in dyy\n',errdyy)    
    fprintf('%5.2e : err in dxxx\n',errdxxx)    
    fprintf('%5.2e : err in dxxy\n',errdxxy)    
    fprintf('%5.2e : err in dxyy\n',errdxyy)    
    fprintf('%5.2e : err in dyyy\n',errdyyy)    
    fprintf('%5.2e : err in dxxxx\n',errdxxxx)    
    fprintf('%5.2e : err in dxxxy\n',errdxxxy)    
    fprintf('%5.2e : err in dxxyy\n',errdxxyy)    
    fprintf('%5.2e : err in dxyyy\n',errdxyyy)    
    fprintf('%5.2e : err in dyyyy\n',errdyyyy)    
end

%%
% some high precision calcs for comparison

% from mathematica
srct = [0;0]; targt = 10^(-4)*[1;1];
kt = sqrt(2);
valt = -0.03670784159258519723+0.24999999750000000625 *1i;
gradxt = -0.00014535819351409084213-0.00002499999987500000021 *1i; 
hessxyt = 0.079577479012801035120+1.249999995833*1e-9*1i;
der3xxyt = 0.000070689659841738687819+0.000012499999916666666823*1i; 
der3xxxt = 1591.54965094568054413183407228+0.00003749999983333333359375*1i;
der4xxyyt = 7.95774786149136280529936105345*1e6+0.12499999875000000364583*1i;
[val,grad,hess,der3,der4] = chnk.helm2d.helmdiffgreen(kt,srct,targt);

abs(val-valt)/abs(valt)
abs(grad(1)-gradxt)/abs(gradxt)
abs(hess(2)-hessxyt)/abs(hessxyt)
abs(der3(1)-der3xxxt)/abs(der3xxxt)
abs(der3(2)-der3xxyt)/abs(der3xxyt) % this derivative appears to be just hard to get
abs(der4(3)-der4xxyyt)/abs(der4xxyyt)
