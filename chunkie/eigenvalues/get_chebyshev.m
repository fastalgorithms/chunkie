function [pol] = get_chebyshev(g, a, b, order)
% get_chebyshev gives the Chebyshev interpolant polynomial of order "order" for
% function "f" on the interval [a,b]

f = @(x) g( (b-a)*x/2 + (b+a)/2 );     % f is interpolant on [-1,1]
n = order;
t = 0:1:n;
x = cos(2*pi*t/(2*n + 1));     
g = f(x);                      
A = cos( (2*pi/(2*n+1))*[0:1:n].'*[0:1:n] );
d = 1/(2*n+1) * A * (g.*[1 2*ones(1,n)]).';
pol = @(t) (cos(t*[0:1:n]).*[1 2*ones(1,n)])*d;
end