function [x,k,g] = gmresF2K(z, K, b, x, opts)
    % Solve (z*I + K) x = b
    % z is a complex number
    % K is a matrix or function handle 
    m = opts.max_iterations;  
    tol = opts.threshold;
    n = numel(b);
    if isa(K,'function_handle')
        Kop = @(v) K(v);
    else
        Kop = @(v) K*v;
    end

    r = b - z*x - Kop(x);
    beta = norm(r); 
    bnorm = norm(b); 
    if bnorm==0, bnorm=1; end
    if beta/bnorm <= tol, return; end
    
    Q = zeros(n,m+1);          % Krylov basis
    H = complex(zeros(m+1,m)); % Hessenberg
    c = complex(zeros(m,1)); 
    s = complex(zeros(m,1));
    g = complex(zeros(m+1,1));
    
    Q(:,1) = r/beta; g(1) = beta;
    kend = m;
    for k = 1:m
        % Arnoldi / MGS
        w = Kop(Q(:,k));
        for i = 1:k
            H(i,k) = Q(:,i)'*w;
            w = w - H(i,k)*Q(:,i);
        end
        H(k+1,k) = norm(w);
        if H(k+1,k) ~= 0
            Q(:,k+1) = w/H(k+1,k); 
        else
            Q(:,k+1) = Q(:,k); 
        end
        
        H(k,k) = H(k,k)+z;

        % Apply previous rotations
        for i = 1:k-1
            Hik  = H(i,k); 
            Hik1 = H(i+1,k);
            H(i,k)   = conj(c(i))*Hik - s(i)*Hik1;
            H(i+1,k) = conj(s(i))*Hik + c(i)*Hik1;
        end

        % New rotation to zero H(k+1,k)
        [c(k),s(k),rlen] = cgivens(H(k,k), H(k+1,k));
        H(k,k) = rlen; H(k+1,k)=0;
    
        % Update g (||r_k|| = |g(k+1)|)
        gk = g(k);
        g(k)   = conj(c(k))*gk - s(k)*g(k+1);
        g(k+1) = conj(s(k))*gk + c(k)*g(k+1);

        % Stop if reached tolerance
        if abs(g(k+1))/bnorm <= tol, kend = k; break; end
    end

    if kend == m
        kend = m;
    end
    
    % Solve using least squares
    y = zeros(kend,1);
    for i = kend:-1:1
        ssum = g(i);
        for j = i+1:kend
            ssum = ssum - H(i,j)*y(j);
        end
        y(i) = ssum / H(i,i);
    end
    
    x = x + Q(:,1:kend)*y;
end

function [c,s,r] = cgivens(a,b)
    % Complex Givens
    % [conj(c) -s; conj(s) c]*[a;b] = [r;0]
    aa = abs(a); bb = abs(b);
    r = hypot(aa,bb);
    if r==0, c=1; s=0; else, c=a/r; s=b/r; end
end