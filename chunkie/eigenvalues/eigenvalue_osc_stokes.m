%%%%%%%%%% Stokes Eigenvalues   %%%%%%%
rad = 1; ctr = [0.0;0.0];
circfun = @(t) ctr + rad*[cos(t(:).');sin(t(:).')];

l = 2.1;
k = 1/(4*l);
h = 10*l; 

rmax = 2/h;
nch = max(10, ceil(2*pi/rmax) + 2);

%%
opts = [];
opts.maxchunklen = 1.1;
pref = []; pref.nchmax = 20000;

chnkr1 = chunkerfuncuni(circfun, nch, opts, pref);
plot(chnkr1, 'k.');
chnkreps = chnkr1;


%%%%%%%% instead of 5 choose any number >=5 for holes in H5 of Vanninathan  

% chnkrs(1,l^2) = chunker();
chunkerhole0 = k*chunkerfuncuni(circfun, 4);
chunkerhole0 = chunkerhole0.reverse();
% nchnkrs = 0;
% nch0 = chunkerhole0.nch;
tic,
nholes = 0;
% nch = chnkreps.nch;
for i = 0:l-1
    for j = 0:l-1
        check1 = [i/l,j/l];
        l1 = norm(check1 + [1/h,1/h]); %% l1 = norm(check1); when only the holes appearing inside those epsilon-cubes which are completely submerged inside Omega is considered 
        l2 = norm(check1 + [1/l,0.0] + [-1/h,1/h]);%% l2 = norm(check1 + [1/l,0.0]); when only the holes appearing inside those epsilon-cubes which are completely submerged inside Omega is considered
        l3 = norm(check1 + [0.0,1/l] + [1/h,-1/h]);%% l2 = norm(check1 + [0.0,1/l]); when only the holes appearing inside those epsilon-cubes which are completely submerged inside Omega is considered
        l4 = norm(check1 + [1/l,1/l] + [-1/h,-1/h]);%% l2 = norm(check1 + [1/l,1/l]); when only the holes appearing inside those epsilon-cubes which are completely submerged inside Omega is considered
        if (l1<=1) && (l2<=1) && (l3<=1) && (l4<=1)
            
            nholes = nholes + 4;
            chnkrctr1 = check1 + (1/(2*l))*[1.0,1.0];
            chnkrctr2 = [-1,1].*chnkrctr1;
            chnkrctr3 = [1,-1].*chnkrctr1;
            chnkrctr4 = [-1,-1].*chnkrctr1;

            chnkrhole1 = chnkrctr1 + chunkerhole0;
            chnkrhole2 = chnkrctr2 + chunkerhole0;
            chnkrhole3 = chnkrctr3 + chunkerhole0;
            chnkrhole4 = chnkrctr4 + chunkerhole0;

            chnkreps = merge([chnkreps , chnkrhole1 , chnkrhole2 , chnkrhole3 , chnkrhole4]);
           
        end
    end
end
%chnkreps = merge([chnkreps, chnkrs(1:nchnkrs)]);
toc
%%
zk = 3.3;
tic, [det, A] = get_determinant(zk, chnkreps); toc
%%

% plot curve, nodes, and normals
plot(chnkr1, 'k.') 
figure(1)
clf
plot(chnkr1,'k.')
hold on
quiver(chnkr1,'k.')
axis equal tight
%%

start = tic;
ncheb = 32;
xx = chnkreps.r(1,:);
yy = chnkreps.r(2,:);
u = (xx.^2 - yy.^2).*sin(xx) + xx.*cos(yy);
v = (xx - yy.^3).*exp(-xx.^2) + yy.*exp(-yy.^2);
% v = @(x) (x(:,1) - x(:,2)^.3).*exp(x(:,1)) + x(:,2).*exp(x(:,1));
             
f1 = chebfun(@(zk) get_determinant(zk, chnkreps), [3,4], ncheb);
[~,A] = chebfun(@(zk) get_determinant(zk, chnkreps), [3,4], ncheb);

%%
plot(real(f1), 'k.');
rts = roots(f1, 'complex');
time1 = toc(start);

% save([num2str(nholes) '.mat'], 'rts', 'time1');

%%

%  g = chebfun(@(x) -besselj(1,x), [0,10]); 
% %%% for Neumann we have to consider the derivative of Bessel function. 
%  rts2 = roots(g, 'complex');


function [det1, varargout] = get_determinant(zk, chnkr1)

    Dk = 2*kernel('ostok', 'd', zk);  
    % for Neumann Dk = 2*kernel('ostok', 'd', zk);...
    % i.e. we have to change sign form - to + in the start of this kernel
    A = chunkermat(chnkr1, Dk);   
    % Builds matrix out of all the points on the boundary, i.e. divide them into panels 

%%
    A = A + eye(2*chnkr1.npt) + normonesmat(chnkr1)/sum(chnkr1.wts(:));
%%

    det1 = det(A);
    varargout{1} = A;
end