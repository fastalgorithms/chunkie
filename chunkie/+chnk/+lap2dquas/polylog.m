function L = polylog(n,z)
%CHNK.LAP2DQUAS.POLYLOG polylogarithm Li_n(z) for positive integer orders
% and complex double z (principal branch, cut along [1,inf)).
%
%   Li_n(z) = sum_{k=1}^inf z^k / k^n,   |z| < 1,
%
% and its analytic continuation elsewhere.
%
% Syntax: L = chnk.lap2dquas.polylog(n,z)
%
% Input:
%   n - array of positive integer orders
%   z - scalar or array of complex numbers (same size as n, or one of the
%       two scalar)
%
% Output:
%   L - Li_n(z), size of the broadcast inputs
%
% Method: direct series for |z| <= 1/2, the inversion formula
%   Li_n(z) + (-1)^n Li_n(1/z) = -(2 pi i)^n/n! B_n(1/2 + log(-z)/(2 pi i))
% for |z| >= 2, and otherwise the expansion in mu = log(z)
%   Li_n(e^mu) = mu^(n-1)/(n-1)! (H_{n-1} - log(-mu))
%                + sum_{k ~= n-1} zeta(n-k) mu^k/k!
% which converges for |mu| < 2 pi.

if isscalar(n) && ~isscalar(z)
    n = n*ones(size(z));
elseif isscalar(z) && ~isscalar(n)
    z = z*ones(size(n));
end
if ~isequal(size(n),size(z))
    error('CHNK.LAP2DQUAS.POLYLOG: n and z must be the same size or scalar');
end
if any(n(:) < 1 | n(:) ~= round(n(:)))
    error('CHNK.LAP2DQUAS.POLYLOG: n must be a positive integer');
end

L = zeros(size(z));
for i = 1:numel(z)
    L(i) = polylog1(n(i),z(i));
end

end

function L = polylog1(n,z)

if n == 1
    L = -log(1-z);
    return
end

az = abs(z);
if az <= 0.5
    L = series(n,z);
elseif az >= 2
    % inversion formula, written in terms of w = 2 pi i x so that the
    % Bernoulli numbers enter only through bounded zeta(2j) values:
    % (2 pi i)^n/n! B_n(x) = w^n/n! - pi i w^(n-1)/(n-1)!
    %                        - 2 sum_{j=1}^{n/2} zeta(2j) w^(n-2j)/(n-2j)!
    w = 1i*pi + log(-z);
    S = w^n/factorial(n) - 1i*pi*w^(n-1)/factorial(n-1);
    for j = 1:floor(n/2)
        S = S - 2*zetaint(2*j)*w^(n-2*j)/factorial(n-2*j);
    end
    L = -(-1)^n*series(n,1/z) - S;
else
    mu = log(z);
    if mu == 0
        L = zetaint(n);
        return
    end
    % terms k = 0...n-2 : zeta(n-k) mu^k / k!
    L = 0;
    for k = 0:n-2
        L = L + zetaint(n-k)*mu^k/factorial(k);
    end
    % k = n-1 : logarithmic term
    H = sum(1./(1:n-1));
    L = L + mu^(n-1)/factorial(n-1)*(H - log(-mu));
    % k = n : zeta(0) = -1/2
    L = L - 0.5*mu^n/factorial(n);
    % k = n-1+2j : zeta(1-2j) = (-1)^j 2 (2j-1)! zeta(2j) / (2 pi)^(2j)
    rat = abs(mu)/(2*pi);
    jmax = ceil(log(eps)/(2*log(rat))) + 2;
    for j = 1:jmax
        k = n-1+2*j;
        c = (-1)^j*2*zetaint(2*j)*exp(gammaln(2*j) - gammaln(k+1));
        L = L + c*mu^(n-1)*(mu/(2*pi))^(2*j);
    end
end

end

function L = series(n,z)
% direct series, |z| <= 1/2
if z == 0
    L = 0;
    return
end
K = ceil(log(eps)/log(abs(z))) + 1;
k = 1:K;
L = sum(z.^k./k.^n);
end

function zv = zetaint(s)
% Riemann zeta at an integer s >= 2
if s >= 30
    k = 1:5;
    zv = sum(k.^(-s));
else
    % psi(m,1) = (-1)^(m+1) m! zeta(m+1)
    zv = (-1)^s*psi(s-1,1)/factorial(s-1);
end
end
