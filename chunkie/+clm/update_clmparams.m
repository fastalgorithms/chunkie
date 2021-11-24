function [clmparams_out] = update_clmparams(clmparams,opts)
    clmparams_out = clmparams;
    ndomain = clmparams.ndomain;
    if(isfield(opts,'rn'))
        clmparams_out.rn = opts.rn;
        
    end
    if(isfield(opts,'lambda'))
        clmparams_out.lambda = opts.lambda;
    end
    rn = clmparams_out.rn;
    lambda = clmparams_out.lambda;
    
    if(strcmpi(clmparams.mode,'te'))
        coef = ones(clmparams.ndomain,1);
    elseif(strcmpi(clmparams.mode,'tm'))
        coef = rn.^2;
    else
        fprintf('Invalid mode specification, reverting to te')
        coef = ones(geom_class.ndomain,1);
    end
    clmparams_out.coef = coef;
    
    
    k = rn/lambda*2*pi;
    kinf = inf*ones(1,ndomain);
    for i=1:ndomain
        if(clmparams.is_inf(i) ~= 0)
            kinf(i) = k(i); 
        end
    end
    c1 = log(1.0/eps)/min(kinf);
    c2 = c1/2;
    a = clmparams.verts(1,clmparams.lvert);
    b = clmparams.verts(1,clmparams.rvert);
    
    xmin = clmparams.xylim(1);
    xmax = clmparams.xylim(2);
    d1 = a-xmin;
    d2 = xmax-b;
    d3 = 2*(b-a);
    C = 2*max([d1 d2 d3]);
    L(1) = C-a;
    L(2) = b+C;
    clmparams_out.cparams{1}.ta = -L(1);
    clmparams_out.cparams{2}.tb = L(2);
    clmparams_out.cpars{1}.L = L(1);
    clmparams_out.cpars{1}.c1 = c1;
    clmparams_out.cpars{1}.c2 = c2;
    clmparams_out.cpars{2}.L = L(2);
    clmparams_out.cpars{2}.c1 = c1;
    clmparams_out.cpars{2}.c2 = c2;
    
    
    k1 = zeros(1,clmparams.ncurve);
    k2 = zeros(1,clmparams.ncurve);
    for i=1:clmparams.ncurve
        k1(i) = k(clmparams.c(1,i));
        k2(i) = k(clmparams.c(2,i));
    end
    clmparams_out.k = k;
    clmparams_out.k1 = k1;
    clmparams_out.k2 = k2;
    
    
    % Determine nch
    nch = zeros(1,clmparams.ncurve);
   

    n0 = 12;
    fac = 1.2/2/pi;
    ppw = 5;

    for i=1:clmparams.ncurve
      ctype = clmparams.cparams{i}.curvetype;
      if(ctype == 0)
          L = clmparams_out.cparams{i}.tb - clmparams_out.cparams{i}.ta;
      elseif(ctype == 1)
          L = sqrt(sum((clmparams_out.cpars{i}.v0-clmparams_out.cpars{i}.v1).^2,1));
      elseif(ctype == 2)
          L = sqrt(sum((clmparams_out.cpars{i}.v1-clmparams_out.cpars{i}.v0).^2,1))/2/sin(clmparams_out.cpars{i}.theta/2);
      elseif(ctype == 3)
          a = clmparams_out.cpars{i}.a;
          b = clmparams_out.cpars{i}.b;
          n = clmparams_out.cpars{i}.n;
          A = clmparams_out.cpars{i}.A;
          alpha = (A*n*pi/(b-a))^2;
          f = @(x) sqrt(1.0 + alpha*cos(x).^2)*(b-a)*n/pi/2;
          L = integral(f,0,pi);
          
      end
      kmax= max(abs(k1(i)),abs(k2(i)));

      nch(i) = round(ppw*L*kmax/32/pi) + n0;
    end
    clmparams_out.nch = nch;
    clmparams_out.npts = sum(nch)*clmparams_out.ngl;
    clmparams_out.nsys = 2*sum(nch)*clmparams_out.ngl;
    
    % diagonal constant for each curve
    alpha1 = zeros(1,clmparams_out.ncurve);
    alpha2 = zeros(1,clmparams_out.ncurve);
    for i=1:clmparams_out.ncurve
      alpha1(i) = 2/(coef(clmparams.c(1,i))+coef(clmparams.c(2,i)));
      alpha2(i) = 2/(1/coef(clmparams.c(1,i))+1/coef(clmparams.c(2,i)));
    end
    clmparams_out.alpha1 = alpha1;
    clmparams_out.alpha2 = alpha2;

end