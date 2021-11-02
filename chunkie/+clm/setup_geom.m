function [clmparams] = setup_geom(geom_class)
    clmparams= [];
    clmparams.xylim = geom_class.xylim;
    clmparams.ndomain = geom_class.ndomain;
    clmparams.rn = geom_class.rn;
    clmparams.lambda = geom_class.lambda;
    ndomain = geom_class.ndomain;
    rn = clmparams.rn;
    lambda = clmparams.lambda;
    
    if(isfield(geom_class,'mode'))
        clmparams.mode = geom_class.mode;
    else
        clmparams.mode = 'te';
    end
    clmparams.verts = geom_class.verts;
    verts = geom_class.verts;
    
    if(strcmpi(clmparams.mode,'te'))
        coef = ones(geom_class.ndomain,1);
    elseif(strcmpi(clmparams.mode,'tm'))
        coef = rn.^2;
    else
        fprintf('Invalid mode specification, reverting to te')
        coef = ones(geom_class.ndomain,1);
    end
    clmparams.coef = coef;
    
    % Determine the extent to go out on either side
    a = verts(1,geom_class.lvert);
    b = verts(1,geom_class.rvert);
    
    xmin = clmparams.xylim(1);
    xmax = clmparams.xylim(2);
    
    is_inf = zeros(1,ndomain);
    k = rn/lambda;
    kinf = inf*ones(1,ndomain);
    for i=1:ndomain
        if(geom_class.region{i}.is_inf ~= 0)
            is_inf(i) = geom_class.region{i}.is_inf;
            kinf(i) = k(i); 
        end
    end
    c1 = log(1.0/eps)/min(kinf);
    c2 = c1/2;
    
    d1 = a-xmin;
    d2 = xmax-b;
    d3 = 2*(b-a);
    C = max([d1 d2 d3]);
    L(1) = C-a;
    L(2) = b+C;
    
    
    
    
    % Determine the number of curves actually used
    % There are two additional curves corresponding to
    % left and right intervals and moreover additional 
    % curves get created if there is a polygon with
    % more than one edge each edge needs to be a separate
    % curve
    
    ncurve = geom_class.ncurve;
    ncurve_use = 2;
    icurve_list = cell(1,ncurve);
    for i=1:ncurve    
        if(geom_class.curves{i}.curvetype == 1)
            iadd = length(geom_class.curves{i}.vert_list)-1;
            iarr = 1:1:iadd;
            icurve_list{i} = ncurve_use + iarr;
        else
            iadd = 1;
            icurve_list{i} = ncurve_use + 1;
        end
        
        ncurve_use = ncurve_use + iadd;
        
    end
    %curve_list = cell(1,ncurve_use);
    cparams = cell(1,ncurve_use);
    cpars = cell(1,ncurve_use);
    
    
    
    % Populate all curve info
    cparams{1}.ta = -L(1);
    cparams{1}.tb = a;
    cparams{1}.curvetype = 0;
    cparams{1}.ilr = 1;
    cparams{1}.ifclosed = false;
    cparams{2}.ta = b;
    cparams{2}.tb = L(2);
    cparams{2}.curvetype = 0;
    cparams{2}.ilr = 2;
    cparams{2}.ifclosed = false;
    cpars{1}.L = L(1); cpars{1}.c1 = c1; cpars{1}.c2 = c2;
    cpars{1}.lvert = geom_class.lvert;
    cpars{2}.L = L(2); cpars{2}.c1 = c1; cpars{2}.c2 = c2;
    cpars{2}.rvert = geom_class.rvert;
    
    
    icurve_cur = 3;
    for i=1:ncurve
        ctype = geom_class.curves{i}.curvetype;
        if(ctype == 1)
            for j=1:length(geom_class.curves{i}.vert_list)-1
                cparams{icurve_cur}.ta = 0;
                cparams{icurve_cur}.tb = 1;
                cparams{icurve_cur}.curvetype = ctype;
                cparams{icurve_cur}.ifclosed = false;
                cpars{icurve_cur}.v0 = verts(:,geom_class.curves{i}.vert_list(j));
                cpars{icurve_cur}.v0ind = geom_class.curves{i}.vert_list(j);
                cpars{icurve_cur}.v1 = verts(:,geom_class.curves{i}.vert_list(j+1));
                cpars{icurve_cur}.v1ind = geom_class.curves{i}.vert_list(j+1);
                icurve_cur = icurve_cur + 1;
            end
        elseif(ctype == 2)
            thet = geom_class.curves{i}.theta;
            cparams{icurve_cur}.ta = -thet/2;
            cparams{icurve_cur}.tb = thet/2;
            cparams{icurve_cur}.ifclosed = false;
            cparams{icurve_cur}.curvetype = ctype;
            cpars{icurve_cur}.v0 = verts(:,geom_class.curves{i}.vert_list(1));
            cpars{icurve_cur}.v0ind = geom_class.curves{i}.vert_list(1);
            cpars{icurve_cur}.v1 = verts(:,geom_class.curves{i}.vert_list(2));
            cpars{icurve_cur}.v1ind = geom_class.curves{i}.vert_list(2);
            cpars{icurve_cur}.theta = geom_class.curves{i}.theta;
            cpars{icurve_cur}.ifconvex = geom_class.curves{i}.ifconvex;
            icurve_cur = icurve_cur+1;
        elseif(ctype == 3)
            nwig = geom_class.curves{i}.nwiggles;
            amp = geom_class.curves{i}.amplitude;
            cparams{icurve_cur}.ta = 0;
            cparams{icurve_cur}.tb = nwig*pi;
            cparams{icurve_cur}.ifclosed = false;
            cparams{icurve_cur}.curvetype = ctype;
            cpars{icurve_cur}.a = verts(1,geom_class.curves{i}.vert_list(1));
            cpars{icurve_cur}.v0ind = geom_class.curves{i}.vert_list(1);
            cpars{icurve_cur}.b = verts(1,geom_class.curves{i}.vert_list(2));
            cpars{icurve_cur}.v1ind = geom_class.curves{i}.vert_list(2);
            cpars{icurve_cur}.A = amp;
            cpars{icurve_cur}.n = nwig;
            icurve_cur = icurve_cur+1;
        end            
    end
    clmparams.cpars = cpars;
    clmparams.cparams = cparams;
    clmparams.ncurve = ncurve_use;
    
    
    % Assemble new clist
    clist = cell(1,ndomain);
    for i=1:ndomain
        for j=1:length(geom_class.region{i}.icurve_list)
            if(geom_class.region{i}.icurve_list(j) < 0)
                clist{i} = [clist{i} -icurve_list{abs(geom_class.region{i}.icurve_list(j))}]; 
            else
                clist{i} = [clist{i} icurve_list{abs(geom_class.region{i}.icurve_list(j))}]; 
            end
        end 
        if(geom_class.region{i}.is_inf == 1)
            clist{i} = [1 clist{i} 2];  
        elseif(geom_class.region{i}.is_inf == -1)
            clist{i} = [-2 clist{i} -1];
        end      
    end
    clmparams.clist = clist;
    
    % Compute domain indices for each curve
    c = zeros(2,ncurve_use);
    for i=1:ndomain
        for j=1:length(clist{i})
            icurve = clist{i}(j);
            if(icurve<0)
                c(2,-icurve) = i;
            else
                c(1,icurve) = i;
            end
        end 
    end
    
    
    k1 = zeros(1,ncurve_use);
    k2 = zeros(1,ncurve_use);
    for i=1:ncurve_use
        k1(i) = k(c(1,i));
        k2(i) = k(c(2,i));
    end
    
    clmparams.k = k;
    clmparams.c = c;
    clmparams.k1 = k1;
    clmparams.k2 = k2;
    clmparams.is_inf = is_inf;
    
    % Get corner info
    [~,ncorner] = size(verts);
    corners = cell(1,ncorner);
    clist = cell(1,ncorner);
    isstart = cell(1,ncorner);
    for i=1:ncorner
        clist{i} = [];
        isstart{i} = [];
    end
    for i=1:ncurve_use
        ctype = cparams{i}.curvetype;
        if(ctype == 1 || ctype == 3)
            v0ind = cpars{i}.v0ind;
            clist{v0ind} = [clist{v0ind} i];
            isstart{v0ind} = [isstart{v0ind} 1];
            v1ind = cpars{i}.v1ind;
            clist{v1ind} = [clist{v1ind} i];
            isstart{v1ind} = [isstart{v1ind} 0]; 
        elseif(ctype == 2)
            ifconvex = cpars{i}.ifconvex;
            if(ifconvex)
                v0ind = cpars{i}.v0ind;
                clist{v0ind} = [clist{v0ind} i];
                isstart{v0ind} = [isstart{v0ind} 1];
                v1ind = cpars{i}.v1ind;
                clist{v1ind} = [clist{v1ind} i];
                isstart{v1ind} = [isstart{v1ind} 0]; 
            else
                v0ind = cpars{i}.v0ind;
                clist{v0ind} = [clist{v0ind} i];
                isstart{v0ind} = [isstart{v0ind} 0];
                v1ind = cpars{i}.v1ind;
                clist{v1ind} = [clist{v1ind} i];
                isstart{v1ind} = [isstart{v1ind} 1];
            end
                
          
        elseif(ctype == 0)
            if(cparams{i}.ilr == 1)
                v0ind = cpars{i}.lvert;
                clist{v0ind} = [clist{v0ind} i];
                isstart{v0ind} = [isstart{v0ind} 0];
            elseif(cparams{i}.ilr == 2)
                v1ind = cpars{i}.rvert;
                clist{v1ind} = [clist{v1ind} i];
                isstart{v1ind} = [isstart{v1ind} 1];
            end
        end 
    end
    for i=1:ncorner
        corners{i}.clist = clist{i};
        corners{i}.isstart = isstart{i};
        corners{i}.nedge = length(clist{i});
    end
    clmparams.issymmetric = 1;
    clmparams.corners = corners;
    clmparams.ncorner = ncorner;
    
    % Determine nch
    nch = zeros(1,ncurve_use);
    fcurve = cell(1,ncurve_use);

    n0 = 12;
    fac = 1.2/2/pi;

    for i=1:ncurve_use
      ctype = clmparams.cparams{i}.curvetype;
      if(ctype == 0)
          L = cparams{i}.tb - cparams{i}.ta;
          fcurve{i} = @(t) clm.complexx6(t,cparams{i}.ilr,cpars{i});
      elseif(ctype == 1)
          L = sqrt(sum((cpars{i}.v0-cpars{i}.v1).^2,1));
          fcurve{i} = @(t) clm.linesegment(t,cpars{i});
      elseif(ctype == 2)
          L = sqrt(sum((cpars{i}.v1-cpars{i}.v0).^2,1))/2/sin(cpars{i}.theta/2);
          fcurve{i} = @(t) clm.circulararc(t,cpars{i});
      elseif(ctype == 3)
          a = cpars{i}.a;
          b = cpars{i}.b;
          n = cpars{i}.n;
          A = cpars{i}.A;
          alpha = (A*n*pi/(b-a))^2;
          f = @(x) sqrt(1.0 + alpha*cos(x).^2)*(b-a)*n/pi/2;
          L = integral(f,0,pi);
          fcurve{i} = @(t) clm.sinecurve(t,cpars{i});
      end
      
      lambda = 1/max(abs(k1(i)),abs(k2(i)));

      nch(i) = round(fac*L/lambda) + n0;
    end
    clmparams.fcurve = fcurve;
    clmparams.nch = nch;
    ngl = 16;
    clmparams.ngl = ngl;
    
    
    clmparams.npts = sum(nch)*ngl;
    clmparams.nsys = sum(nch)*ngl*2;
    
    
    
% diagonal constant for each curve
    alpha1 = zeros(1,ncurve_use);
    alpha2 = zeros(1,ncurve_use);
    for i=1:ncurve_use
      alpha1(i) = 2/(coef(c(1,i))+coef(c(2,i)));
      alpha2(i) = 2/(1/coef(c(1,i))+1/coef(c(2,i)));
    end
    clmparams.alpha1 = alpha1;
    clmparams.alpha2 = alpha2;

    % Set src_in if it exists
    src_in = zeros(2,ndomain);
    issrc_in = false(1,ndomain);
    for i=1:ndomain
        if(isfield(geom_class.region{i},'src_in'))
            src_in(:,i) = geom_class.region{i}.src_in(:);
            issrc_in(i) = true;
        end
    end
    clmparams.src_in = src_in;
    clmparams.issrc_in = issrc_in;
    
    
    
    
end