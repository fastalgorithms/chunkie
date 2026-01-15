%%%%%%%%%%%%%%%%% DESCRIPTION OF DOMAIN %%%%%%%%%%%%%%%%%%%%
%%%%%circfun = @(t) ctr + rad*[cos(t(:).');sin(t(:).')];
%%%%%%%%%%    PERIODIC HOLES            %%%%%%%%%%%%%
%%%%%%%%%%    CODE IS RUNNING FOR N<8   %%%%%%%%%%%%%%
%%%%%%%%%%    CONSIDERING EACH HOLE IS COVERED BY SUB-EPSILON CUBE, 
%%%%%%%%%%    WHICH IS INSIDE AN EPSILON-CUBE, SUCH THAT 
%%%%%%%%%%    SUB-EPSILON CUBE DOESN'T INTERSECT BOUNDARY %%%%%%%%%%    

rad = 1; ctr = [0.0;0.0];
circfun = @(t) ctr + rad*[cos(t(:).');sin(t(:).')];

l = 3;
k = 1/(4*l);
h = 5*l; 

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

tic,
[eigval, ddets, ccoef] =  obj_fun_flam(chnkreps, 10, 15, 24, opts);
toc
%%
save([num2str(nholes) '.mat'], 'eigval', 'ddets', 'ccoef', 'chnkreps');



%%%%%%%%%%%%%%%% FAST CODES %%%%%%%%%%%%%%%%%%%%%%%

function [eigval,varargout] = obj_fun_flam(chnkr0, amin, bmin, ncheb, opts)
                                        
%
%  This function computes the laplacian eigenvalue on the interval
%  [amin, bmin] using ncheb points and returns, 
%  the value of the objective function val, the eigenvalue zk, 
%  dudn given by sigma, mu which is the null vector of -2*D, bie_norm
%  which is the normalization factor, and F which is the flam compressed
%  representation at the eigenvalue.
%
%  The origin is at the point tn \in [-1,1] on chunk ichn
%
%  The routine can also be operated to return val, sig, mu, bie_norm, F
%  if the eigenvalue is known.
%
%

    
    eps = 1e-7;
    if isfield(opts, 'eps')
        eps = opts.eps;
    end
    
    opts_flam = [];
    opts_flam.eps = eps;
    opts_flam.flamtype = 'rskelf';
    opts_flam.forceproxy = true;
    opts_flam.occ = 200;
    
    dval = 1.0;
    if ~isfield(opts, 'zk')
        x0 = cos(pi*(0:(ncheb-1))/(ncheb-1));
        xcheb = (bmin-amin)*(1+x0)/2+amin;
        ddets = zeros(size(xcheb));
        for ii=1:ncheb
            zk = xcheb(ii);
            Sp = 2*kernel('helmholtz','d',zk); 
            % ORIGINALLY THE TYPE OF KERNAL WAS 'sprime' BUT I WILL BE 
            % CHOOSING 'd' DOUBLE-LAYER AND FOR NEUMANN SIGN WILL BE CHANGED TO + FROM -
            F  = chunkerflam(chnkr0, Sp, dval, opts_flam);
            ddets(ii) = exp(rskelf_logdet(F));
        end
    
        ddets = ddets/max(abs(ddets));
        varargout{1} = ddets;
        
    
        %%
        ns = 0:(ncheb-1);
        [NT,XT] = meshgrid(ns,x0);
        tmat = cos(NT.*acos(XT));
        ccoef = tmat\ddets.';
        varargout{2} = ccoef;
        
        
        %%
        cmat = spdiags(ones(2,ncheb).'/2,[-1,1],ncheb,ncheb);
        cmat(1,2) = 1/sqrt(2);
        cmat(2,1) = 1/sqrt(2);
        vvec = ccoef;
        vvec(1) = sqrt(2)*vvec(1);
        cmat(end,:) = cmat(end,:)-1/2*vvec.'/vvec(end);
        es = eigs(cmat,ncheb);
        es(abs(real(es))>0.9) = [];
        es(abs(imag(es))>1E-3) = [];
        eigval = real(es);
        eigval = (bmin-amin)*(eigval+1)/2+amin;
        
        err_eig_real = max(abs(imag(es)));
        varargout{3} = err_eig_real;
        
       
    else
        zk = opts.zk;
    end
end