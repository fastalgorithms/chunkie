% maybe m is not needed
m = 2;
chnkr = get_chunker(1.5);
xx = chnkr.r(1,:);
yy = chnkr.r(2,:);
u = (xx.^2 - yy.^2 + yy).*sin(xx) + cos(xx).*(0.5*yy.^2 + 0.3*yy);
v = (xx - yy.^3 + xx.^2).*cos(-yy) + sin(yy).*(xx.^2 + xx);
A = @(zk) get_matrix(zk, chnkr);
g = @(zk) ones(1, length(u))*( ( (u.').*( A(zk)\(v.') ) ).*(chnkr.wts(:)) );
ff = chebfun(@(zk) 1/g(zk), [3,4]);
% A = @(zk) get_matrix(zk, chnkr);
% g = @(zk) ones(1, length(u))*( ( (u.').*( get_matrix(zk, chnkr)\(v.') ) ).*(chnkr.wts(:)) );
% ff = chebfun(@(zk) 1/g(zk));
ff(3)
return

% f1 = @(j,zk) get_determinant(zk, get_chunker(j)); %% initial complex
% muller argument
% x = zeros(m,1);
% y = zeros(m,1); 
% niters = zeros(m,1);
% 

% f1 = @(j,zk) get_determinant(zk, get_chunker(j));
x = zeros(m,1);
y = zeros(m,1); 
% niters = zeros(m,1);

% for j = 1:m
    %chnkr = get_chunker(j);
    % f = @(zk) get_determinant(zk, chnkr);
    x0 = 3.23537; 
    x1 = 3.23637;
    x2 = 3.23527;

    % if j==1
    %      x0 = 3.23537; 
    %      x1 = 3.24537;
    %      x2 = 3.23837;
    % else
    %     x0 = x(j-1) - 0.02;
    %     x1 = x(j-1);
    %     x2 = x(j-1) + 0.02;
    % end
    y0 = ff(x0); y1 = ff(x1); y2 = ff(x2);
    n = 5; % number of iterations
    for i = 1:n
        h0 = x1 - x0; 
        h1 = x2 - x1;
        delta0 = (y1 - y0)/h0;
        delta1 = (y2 - y1)/h1;
        a = (delta1 - delta0)/(h1 + h0); 
        b = a*h1 + delta1;
        c = y2;

        x3 = x2 - 2*c/(b + sign(b)*sqrt(b*b - 4*a*c));

        x0 = x1; y0 = y1;
        x1 = x2; y1 = y2;
        x2 = x3; y2 = ff(x2);
        if (abs(y2) < 1e-12)
            break;
        end
    end
    niters = i;
% niters(j) = i;
% x(j) = x3;
% y(j) = ff(x3);
% end

% plot(imag(ff), 'k.');
% rts = roots(ff, 'complex');
% 
% %%
% Amat = get_matrix(akk, chnkr); 

function [A] = get_matrix(zk, chnkr1)

    Dk = 2*kernel('helm', 'd', zk);  
    A = chunkermat(chnkr1, Dk);   
%%
    A = A + eye(chnkr1.npt);
end
