function [x,iter,resvec,resfin] = gmresF2K(K, b, z, x, opts)
%GMRESF2K   Stall-free Generalized Minimum Residual Method 
%           for operators of the form A*x = (z*I+K)*x = b,
%           where z is a complex number, and
%           K is a matrix or function handle.
%   returns X: the solution,
%           ITER: number of iterations,
%           RESVEC: residue for all iterations,
%           RESFIN: final value for residue.
%
%   X = GMRES(K,B) attempts to solve the system of linear equations (I+K)*X = B
%   for X.  The N-by-N matrix K must be square and the right hand side 
%   column vector B must have length N. This uses the unrestarted
%   method with MIN(N,100) total iterations.
%
%   X = GMRES(KFUN,B) accepts a function handle KFUN instead of the matrix
%   K. KFUN(X) accepts a vector input X and returns the matrix-vector
%   product K*X. In all of the following syntaxes, you can replace K by
%   KFUN.
%
%   X = GMRES(KFUN,B,Z) attempts to solve the system of linear equations 
%   (Z*I+K)*X = B for X, where Z is a complex number.
% 
%   X = GMRES(KFUN,B,Z,X0) specifies the first initial guess. 
%   If X0 is [] then GMRESF2K uses the default, an all zero vector.
%
%   X = GMRES(KFUN,B,Z,X0,OPTS) specifies the tolerance and maximum 
%   number of iterations of the method.  
%   If OPTS.TOL is [] then GMRESF2K uses the default, 1e-16.
%   If OPTS.MAXIT is [] then GMRESF2K uses the default, min(N,100).
%
%   Example:
%      n = 21; A = gallery('wilk',n);  b = sum(A,2);
%      tol = 1e-12;  maxit = 15; M = diag([10:-1:1 1 1:10]);
%      opts = []; opts.max_iterations = 1e4; opts.threshold = 1e-100;
%      x = gmresF2K(A,b,1,[],opts);
%   Or, use this matrix-vector product function
%      %-----------------------------------------------------------------%
%      function y = kfun(x,n)
%      y = [0; x(1:n-1)] + [((n-1)/2:-1:0)'; (1:(n-1)/2)'].*x+[x(2:n); 0];
%      %-----------------------------------------------------------------%
%   as inputs to GMRESF2K:
%      x1 = gmresF2K(@(x)kfun(x,n),b,1,[],opts);

    if nargin < 2, error('gmresF2K: need at least K and b'); end
    n = numel(b);

    if nargin < 5 || isempty(opts), opts = struct; end
    if ~isfield(opts,'maxit') && ~isfield(opts,'max_iterations')
        opts.maxit = min(n,100);
    end
    if ~isfield(opts,'tol') && ~isfield(opts,'threshold')
        opts.tol = 1e-16;
    end

    if nargin < 4 || isempty(x), x = zeros(size(b)); end
    if nargin < 3 || isempty(z), z = 1; end

    if isa(K,'function_handle')
        Kop = @(v) K(v);
    else
        Kop = @(v) K*v;
    end

    if isreal(z)
        usecomplex = false;
    else
        usecomplex = true;
    end

    n = length(b);
    if isempty(x), x = zeros(size(b)); end
    m = opts.maxit;
    tol = opts.tol;

    r = b-z*x-Kop(x);
    
    beta = norm(r);
    b_norm = norm(b);

    % Preallocate
    Q = zeros(n, m+1);
    H = zeros(m+1, m);
    cs = zeros(m,1);
    sn = zeros(m,1);
    g = zeros(m+1,1);
    resvec = zeros(m,1);
    % Init
    Q(:,1) = r / beta;
    g(1) = beta;

    if b_norm < tol || abs(g(1)) / b_norm < tol
        return
    end
    
    for k = 1:m

        % Arnoldi
        [H(1:k+1, k), Q(:, k+1)] = arnoldi(K, Q, k);
        H(k, k) = H(k, k) + z;
        [H(1:k+1, k), cs(k), sn(k)] = apply_givens_rotation(H(1:k+1,k), cs, sn, k, usecomplex);
        
        % Apply Givens rotation 
        g(k + 1) = -sn(k) * g(k);
        g(k)     =  cs(k) * g(k);
        
        % Find residue and error 
        %err(k) = norm(Q(:, 1:k+1)*g(k + 1));
        
        yk = H(1:k, 1:k) \ g(1:k);
        xk = x + Q(:, 1:k) * yk;
        resvec(k) = norm(b-z*xk-Kop(xk)); 
        
        current_relative_tol = abs(g(k + 1)) / b_norm;
        if (current_relative_tol <= tol)
            break;
        end

    end
    
    y = H(1:k, 1:k) \ g(1:k);
    x = x + Q(:, 1:k) * y;
    iter = k;
    resfin = norm(b-z*x-Kop(x));
    resvec = resvec(1:k);
end

function [h, cs_k, sn_k] = apply_givens_rotation(h, cs, sn, k, usecomplex)
    for i = 1:k-1
        temp   =  cs(i)*h(i) + sn(i)*h(i+1);
        h(i+1) = -sn(i)*h(i) + cs(i)*h(i+1);
        h(i)   = temp;
    end
    if usecomplex
        rotate = @cgivens_rotation;
    else
        rotate = @givens_rotation;
    end
    [cs_k, sn_k] = rotate(h(k), h(k+1));
    h(k)   = cs_k*h(k) + sn_k*h(k+1);
    h(k+1) = 0.0;
end

function [cs, sn] = givens_rotation(v1, v2)
    if (v1 == 0)
        cs = 0;
        sn = 1;
    else
        t = sqrt(v1^2 + v2^2);
        cs = v1 / t; 
        sn = v2 / t;
    end
end

function [c,s,r] = cgivens_rotation(a,b)
    % Complex Givens: [conj(c) -s; conj(s) c]*[a;b] = [r;0], 
    aa = abs(a); bb = abs(b);
    r = hypot(aa,bb);
    if r==0, c=1; s=0; else, c=a/r; s=b/r; end
end