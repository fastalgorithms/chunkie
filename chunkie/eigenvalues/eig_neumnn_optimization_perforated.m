%%%%%%%%%%  COMPLEX MULLER FOR ROOTS OF FREDHOLM DETERMINANT %%%%%%%%%

%%%%%%%%%%%% eigval computation for neumann laplacian using the complex muller
% m = 2;

% f1 = @(j,zk) get_determinant(zk, get_chunker(j));
% x = zeros(m,1);
% y = zeros(m,1); 
% niters = zeros(m,1);
% 
% for j = 1:m
%     chnkr = get_chunker(j);
%     f = @(zk) get_determinant(zk, chnkr);
%     x0 = 3.81; 
%     x1 = 3.87;
%     x2 = 3.88;
% 
%     if j==1
%         x0 = 3.81; 
%         x1 = 3.87;
%         x2 = 3.88;
%     else
%         x0 = x(j-1) - 0.02;
%         x1 = x(j-1);
%         x2 = x(j-1) + 0.02;
%     end
%     y0 = f(x0); y1 = f(x1); y2 = f(x2);
%     n = 50; % number of iterations
%     for i = 1:n
%         h0 = x1 - x0; 
%         h1 = x2 - x1;
%         delta0 = (y1 - y0)/h0;
%         delta1 = (y2 - y1)/h1;
%         a = (delta1 - delta0)/(h1 + h0); 
%         b = a*h1 + delta1;
%         c = y2;
% 
%         x3 = x2 - 2*c/(b + sign(b)*sqrt(b*b - 4*a*c));
% 
%         x0 = x1; y0 = y1;
%         x1 = x2; y1 = y2;
%         x2 = x3; y2 = f(x2);
%         if (abs(y2) < 1e-12)
%             break;
%         end
%     end
% niters(j) = i;
% x(j) = x3;
% y(j) = f(x3);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolating the optimization function requires around 300 degree
% polynomials
% ncheb = 30;
chnkr = get_chunker(1.5);
xx = chnkr.r(1,:);
yy = chnkr.r(2,:);
u = (xx.^2 - yy.^2 + yy).*sin(xx) + cos(xx).*(0.5*yy.^2 + 0.3*yy);
v = (xx - yy.^3 + xx.^2).*cos(-yy) + sin(yy).*(xx.^2 + xx);
A = @(zk) get_matrix(zk, chnkr);
g = @(zk) ones(1, length(u))*( ( (u.').*( A(zk)\(v.') ) ).*(chnkr.wts(:)) );
ff = chebfun(@(zk) 1/g(zk), [3.4,3.5]);
%%
plot(imag(ff), 'k.');
rts = roots(ff, 'complex');

%%
%Amat = get_matrix(akk, chnkr); 

function [A] = get_matrix(zk, chnkr1)

    Dk = 2*kernel('helm', 'd', zk);  
    A = chunkermat(chnkr1, Dk);   
%%
    A = A + eye(chnkr1.npt);
end


%%%%%%%%% OBSERVATIONS FROM THIS CODE %%%%%%%%
% 1. For the number of holes = 12 and for the function 
%    ff = chebfun(@(zk) 1/g(zk), [3.4,3.5]); 
%    we see that only 18 degree polynomials is needed,
%    which means that if we have made precise the interval where 
%    we are looking for the eigenvalues then  we are able to use
%    this function as the optimization function and try it
%    to calculate the eigenvalue of the operator.
%    Also, ur agrees with the Fredholm determinant technique
% 2. For the number of holes = 16 and for the function 
%    ff = chebfun(@(zk) 1/g(zk), [3.4,3.5]);
%    we see that only 18 degree polynomials is needed and the roots agrees
%    with the Fredholm techniques.
% 3. For the number of holes = 24 and for the function 
%    ff = chebfun(@(zk) 1/g(zk), [3.4,3.5]);
%    we see that only 13 degree polynomials is needed and the roots agrees
%    with the Fredholm techniques, while for the Fredholm determinant we need 
%    15 degree polynomial. It appears that the optimization function may
%    be needing lesser degree interpolant as compared to the Fredholm as
%    the number of holes grows.