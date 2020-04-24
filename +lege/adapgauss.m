function [rint,maxrec,numint,ier] = adapgauss(fun,a,b,t,w)

k = 16;
eps = 1e-12;

if nargin < 5
    [t,w] = lege.exps(k);
end

nnmax=100000;
maxdepth=200;

stack = zeros(2,maxdepth);
vals = zeros(maxdepth);

% start the recursion

stack(1,1)=a;
stack(2,1)=b;

vals(1) = oneintp(fun,a,b,t,w);

% recursively integrate the thing

j=1;
rint=0;
ier=0;
maxrec=0;
for i=1:nnmax
    numint=i;
    if(j > maxrec); maxrec=j; end

%       subdivide the current subinterval

    c=(stack(1,j)+stack(2,j))/2;
    v2 = oneintp(fun,stack(1,j),c,t,w);
    v3 = oneintp(fun,c,stack(2,j),t,w);
    
    dd= abs(v2+v3-vals(j));
    if(dd <= eps) 


%       if the function on this subinterval has been 
%       integrated with sufficient accuracy - add the 
%       value to that of the global integral and move up
%       in the stack

%
        rint=rint+v2+v3;
        j=j-1;
%
%        if the whole thing has been integrated - return
%
        if(j == 0); return; end

    else

%       if the function on this subinterval has not been 
%       integrated with sufficient accuracy - move 
%       down the stack

        stack(1,j+1)=stack(1,j);
        stack(2,j+1)=(stack(1,j)+stack(2,j))/2;
        vals(j+1)=v2;
%
        stack(1,j)=(stack(1,j)+stack(2,j))/2;
        vals(j)=v3;
    %
        j=j+1;
    %     
    %       if the depth of the recursion has become excessive - bomb
    %
        if(j > maxdepth) 
            ier = 8;
            return;
        end 
    end
end

ier = 16;
end

function val = oneintp(fun,a,b,t,w)
%       integrate the function fun on the interval [a,b]
 
u=(b-a)/2;
v=(b+a)/2;
tt = u*t+v;
ftt = fun(tt);
val = sum(ftt.*w)*u;
end

