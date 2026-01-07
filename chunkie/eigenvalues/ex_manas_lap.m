%%%%%% Disretising the Boundary Integral Equation for the equation 
%%%%%% lap(u) = 0     on the domain B(0,1)
%%%%%% u = f       on the boundary of B(0,1)
clear;

%%%%%%%% Definition of function and the boundary  %%%%%%%%%%

x0 = [5.1 3.14];
f = @(x) 0.5*log( (x(:,1) - x0(:,1)).^2 + (x(:,2) - x0(:,2)).^2 );
a = 0.4;m = 3;b = 2;

% ge = @(t) [a*cos(t);b*sin(t)];
% der_ge = @(x) [-a*sin(x); b*cos(x)];
% der2_ge = @(x) [-a*cos(x); -b*sin(x)];
% nge = @(t) [b*cos(t);a*sin(t)]./( ((b*cos(t)).^2 + (a*sin(t)).^2).^0.5 );

ge = @(t) (1 + a*cos(m*t)).*[cos(t);sin(t)];
der_ge = @(t) -a*m*sin(m*t).*[cos(t);sin(t)] + (1 + a*cos(m*t)).*[-sin(t);cos(t)];
der2_ge = @(t) -(1 + a*cos(m*t)).*[cos(t);sin(t)] + 2*a*m*sin(m*t).*[sin(t);-cos(t)]...
            -a*m*m*cos(m*t).*[cos(t);sin(t)];
speed =  @(t) -a*m*sin(m*t).*[sin(t);-cos(t)] + (1 + a*cos(m*t)).*[cos(t);sin(t)];
nge = @(t) speed(t)./vecnorm(speed(t));


%%%%  Defn of pts where diff betweenn true sol and double_layer is evaluated  %%%%

jj = 100;
pp = -pi:pi/jj:pi;
dir = [cos(pp);sin(pp)].';
rng(1);                                    %%%%%% random number generator fixes seed
rr = 0.5 * min(a,b) * rand(length(pp),1);
x = rr.*dir;
xp = x(:,1); yp = x(:,2);

%%%%% Defn of true sol  %%%%%%%

sol = f(x);    



n = 100;
t = -2*pi*n/(2*n+1) : 2*pi/(2*n+1) : pi;

%%%% Calculating bdry pts and bdry normal for sigma evaluation  %%%%
xx = ge(t).';
nxx = nge(t).';
    
l = length(xx(:,1));
y = f(xx);           % y for f = log|x-x0|   or   f = 5
w = (2*pi)/(2*n+1);
absder_g = (((der_ge(t).').^2)*[1;1]).^(0.5);
ker = get_kerneldl(t, ge, nge, der_ge, der2_ge);
A = -0.5*eye(2*n+1) + w*ker.*(absder_g.');
sigma = A\y;
xg = xx(:,1).'; yg = xx(:,2).';
xn = nxx(:,1).'; yn = nxx(:,2).';
xgm = repmat(xg, length(xp), 1);
ygm = repmat(yg, length(yp), 1);

%%%%%  Calculation of Sigma %%%%%%
ns = 1:1:100;
err = zeros(length(ns), 1);
for ii = 1:length(ns)
    % n = ns(ii);
    % n = 100;
    % t = -2*pi*n/(2*n+1) : 2*pi/(2*n+1) : pi;
    % 
    % %%%% Calculating bdry pts and bdry normal  %%%%
    % xx = ge(t).';
    % nxx = nge(t).';
    % 
    % l = length(xx(:,1));
    % y = f(xx);           % y for f = log|x-x0|   or   f = 5
    % w = (2*pi)/(2*n+1);
    % absder_g = (((der_ge(t).').^2)*[1;1]).^(0.5);
    % ker = get_kerneldl(t, ge, nge, der_ge, der2_ge);
    % A = -0.5*eye(2*n+1) + w*ker.*(absder_g.');
    % sigma = A\y;


    jj = 100;
    pp = -pi:pi/jj:pi;
    dir = [cos(pp);sin(pp)].';
    rng(1);                                    %%%%%% random number generator fixes seed
    rr = (1.5*ii/101) * min(a,b) * rand(length(pp),1);
    x = rr.*dir;
    xp = x(:,1); yp = x(:,2);

    %%%%% Defn of true sol  %%%%%%%

    sol = f(x);    

    %%%%%%% Calculating Double Layer %%%%%%%
    
   
    % xg = xx(:,1).'; yg = xx(:,2).';
    % xn = nxx(:,1).'; yn = nxx(:,2).';
    xpt = repmat(xp, 1, length(xg)); 
    ypt = repmat(yp, 1, length(yg)); 
    double_layer = (w/(2*pi))*( ( (xpt - xgm).*xn + (ypt - ygm).*yn )./( ...
                        (xpt - xgm).^2 + (ypt - ygm).^2 ) )*(sigma.*absder_g);

    
    %%%%%%%%%% Calculating error between true sol and Double Layer %%%%%%%%%
                     
    err(ii) = max(abs(sol - double_layer));         
end
figure; semilogy(1.5*ns/101, err, 'r.'); 
% hold on; semilogy(ns, exp(-0.3*ns), 'b.')




function [ker] = get_kerneldl(t, gamma, normal, der1, der2)
x = gamma(t).';
n = normal(t).';
d1 = der1(t).';
d2 = der2(t).';
xx = x(:,1); xy = x(:,2);
nx = n(:,1).'; ny = n(:,2).';
xpt = repmat(xx, 1, length(xx)); xgm = repmat(xx.', length(xx), 1);
ypt = repmat(xy, 1, length(xy)); ygm = repmat(xy.', length(xy), 1);
ker = (1/(2*pi))*((xpt-xgm).*nx + (ypt-ygm).*ny)./((xpt-xgm).^2 + (ypt-ygm).^2);
for i = 1:length(ker(:,1))
    ker(i,i) = (1/(4*pi))*( n(i,1)*d2(i,1) + n(i,2)*d2(i,2) )/(  d1(i,1)^2 + d1(i,2)^2  );
    % ker(i,i) = -1/(4*pi); % for circle
end
end






% % % % % 
% % % % % 
% % % % % for m = 1:701
% % % % % % m = -56;
% % % % % f = @(t) exp(1j*m*t);
% % % % % % gamma = @(x) [cos(x);sin(x)];
% % % % % % der_gamma = @(x) [-sin(x);cos(x)];
% % % % % % jj = 501;
% % % % % % pp = -pi:pi/jj:pi;
% % % % % % dir = [cos(pp);sin(pp)].';
% % % % % % rng(1);                                    %%%%%% random number generator fixes seed
% % % % % % rr = 0.9*rand(length(pp),1);
% % % % % % x = rr.*dir;
% % % % % % n = 501;
% % % % % % t = -2*pi*n/(2*n+1) : 2*pi/(2*n+1) : pi;
% % % % % % xx = gamma(t).';
% % % % % % l = length(xx(:,1));
% % % % % y = -2*f(t).';
% % % % % % y = -2*f(x);
% % % % % w = 2*pi/(2*n+1);
% % % % % ker = 1/(2*pi);
% % % % % A = eye(2*n+1) + ker*w*ones(2*n+1);
% % % % % sigma = A\y;
% % % % % double_layer = zeros(1,length(pp));
% % % % % for kk = 1:length(pp)
% % % % %    jl = 0;
% % % % %    for kkk = 1:length(xx)
% % % % %        jl = jl + (w/(2*pi)) * sigma(kkk) * ( (x(kk,1) - xx(kkk,1))*xx(kkk,1) + ....
% % % % %                 (x(kk,2) - xx(kkk,2))*xx(kkk,2) )/( (x(kk,1) - xx(kkk,1))^2 + ....
% % % % %                 (x(kk,2) - xx(kkk,2))^2 ) ;
% % % % %    end
% % % % %     double_layer(kk) = jl;
% % % % % end
% % % % % % xp = x(:,1); yp = x(:,2);
% % % % % % xg = xx(:,1).'; yg = xx(:,2).';
% % % % % % xpt = repmat(xp, 1, length(xg)); xgm = repmat(xg, length(xp), 1);
% % % % % % ypt = repmat(yp, 1, length(yg)); ygm = repmat(yg, length(yp), 1);
% % % % % 
% % % % % dott = (w/(2*pi))*( ( (xpt - xgm).*xg + (ypt - ygm).*yg )./( ...
% % % % %     (xpt - xgm).^2 + (ypt - ygm).^2 ) )*sigma;
% % % % % 
% % % % % 
% % % % % 
% % % % % sol = f(pp).*((rr.^abs(m)).');
% % % % % % sol = 5;
% % % % % err(m) = max(abs(sol - double_layer));
% % % % % % err(m) = max(abs(sol.' - dott));
% % % % % % err1 = max(abs(sol - doubl_layer));
% % % % % end
% function [val] = boundary_term(x)
% a = x(1); b = x(2);
% val = a^2;
% end

%%%%%%%%%%%    Learning points
% if your discretising points on "boundary" are all distant 'h' away from each other,
% then you can get the error around points in the "domain" which are at atmost 
% distance 5h from the boundary. 


%%%%%%%%   Homeworks %%%%%%
% 1. Fix an 'n' and Vectorise code 
% 2. consider f = e^{ i m\theta} for m \neq 0, sol = r^(|m|).e^{ i m\theta}  
% 3. consider the following f = -log|x-x_0| where |x_0|>1
% 4. Discretize [x(t) y(t)] = (1+a.cos(mt))[cos(t) sin(t)]

% 5. Plot the error as a function of the number 'n' and try to get the error 
%    and try to justify it in the case of Laplacian.
%   ANSWER:- For a fixed x(according to our notations which is the pts inside domain 
%            where we evaluate the error), the error behaves as the exp(-c.n) where c
%            depends upon distance between the boundary and x, as the kernel
%            is of the order 1/(x-y) when x is in domain and y is in boundary.


% 6. Now try to solve the Helmholtz equation and see the asymptotics of
%    Single Layer and Double Layer in this case, where Green's function is
%    given by G(x,y) = i/4*Hankel(0,1)(k|x-y|), i.e., calculate the limit as x
%    goes to y of G_s(x,y) and limit as x goes to y of G_d(x,y) 
%   ANSWER:- Look at the notes titled NODES exercise section 

% 7. Plot the error in this as done for the Laplacian in the above code and 
%    get the rates and explain the reason for the deacys. 

%%%%%%%   Observations %%%%%%
% 1. for the case of laplacian green's function, error is given by
%    err(n) = exp(-0.68 * n), which we have observed by fixing the pts(x in
%    our notations), where error is to be evaluated, and varying the nodes.
% 2. for the case of laplacian green's function, error depends on the
%    distance of interior points and boundary points, which is the "5h rule"
%    which we have observed by fixing the node as 100, and varying the distance of 
%    pts(x in our notations), where error is to be evaluated, from the boundary 
% 3. The nearer the points(x in our notations), where we are evaluating the difference
%    (u(x) - double_layer(x)), the more the error, which is basically the point number 2 