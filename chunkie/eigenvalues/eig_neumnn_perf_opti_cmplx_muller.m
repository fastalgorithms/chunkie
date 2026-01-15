% f = @(x) sin(x);
% g = @(x) cos(x);




% maybe m is not needed
% m = 2;
chnkr = get_chunker(3);
xx = chnkr.r(1,:);
yy = chnkr.r(2,:);
u = (xx.^2 - yy.^2 + yy).*sin(xx) + cos(xx).*(0.5*yy.^2 + 0.3*yy);
v = (xx - yy.^3 + xx.^2).*cos(-yy) + sin(yy).*(xx.^2 + xx);
A = @(zk) get_matrix(zk, chnkr);
g = @(zk) ones(1, length(u))*( ( (u.').*( A(zk)\(v.') ) ).*(chnkr.wts(:)) );
ff = @(zk) 1/g(zk);

x0 = 3.409274696865717; x1 = 3.409274686865717; x2 = 3.409274686865787;
y0 = ff(x0); y1 = ff(x1); y2 = ff(x2);
n = 10; % number of iterations
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
 

function [A] = get_matrix(zk, chnkr1)

    Dk = 2*kernel('helm', 'd', zk);  
    A = chunkermat(chnkr1, Dk);   
%%
    A = A + eye(chnkr1.npt);
end
