function [clmparams] = setup_geom(geom_class)
%
% clmparams stores the following quantities:
%   * xylim = min and max (x,y) coorindates where solution is desired
%               if manually specified, must completely contain
%               non flat regions of the geometry
%   * ndomain = number of regions in the geometry
%   * rn = refractive index of each region
%   * lambda = wavelength of incident light
%   * mode = Whether to use 'te' or 'tm' mode. 
%       Here we always assume that the transmission problem is either in
%       the 'TE' or 'TM' polarization, which implies that the boundary conditions
%       are of the form
%       [u] = f, [coef du/dn] = g
%
%       where [u] denotes the jump in u across an interface. 
%       To obtain the TE polarization, coef \equiv 1, and for the 
%         TM polarization
%         coef should be set to the square of the refractive index of 
%         the medium.
%    * verts = vertices defining the geometry (these are the corner 
%         locations at which RCIP will be used)
%    * coef (ndomain,1) = inferred from rn, lambda, and mode. Determines
%          coef int he boundary condition above.
%    * lvert,rvert = vertex number corresponding to left/right most vertex beyond
%           which complex coordinates are used
%    * ngr = grid size for plotting solution
%    * cpars,cparams = parameters defining the curves, cell array of structs
%               cparams contains left and right end points of parameter
%               space, whether the curve is closed or not and the type
%               of the curve
%    * ncurve = number of curves
%    * clist = list of size (ndomain,1), ordered list of curves defining
%          the the regions.
%    * k = (ndomain,1) wavenumber in each region and is given by
%         rn*2*pi/lambda
%    * c = (2,ndomain) for each curve speficies which region is in the
%         positive normal direction and negative normal direction
%    * k1,k2 (1,ndomain) = k1 = k(c(1,i)); and k2 = k(c(2,i)) wave number
%          on the positive and negative normal side of each curve
%    * is_inf 

    clmparams= [];
    
    
    if(isfield(geom_class,'xylim'))
        clmparams.xylim = geom_class.xylim;
        
        xmin = geom_class.xylim(1);
        xmax = geom_class.xylim(2);
        ymin = geom_class.xylim(3);
        ymax = geom_class.xylim(4);
        
        xmin1 = min(geom_class.verts(1,:));
        xmax1 = max(geom_class.verts(1,:));
        ymin1 = min(geom_class.verts(2,:));
        ymax1 = max(geom_class.verts(2,:));
        
        
        xmin = min(xmin,xmin1);
        xmax = max(xmax,xmax1);
        ymin = min(ymin,ymin1);
        ymax = max(ymax,ymax1);
        
    else
        xmin = min(geom_class.verts(1,:));
        xmax = max(geom_class.verts(1,:));
        ymin = min(geom_class.verts(2,:));
        ymax = max(geom_class.verts(2,:));
    end
    
    dx = xmax - xmin;
    dy = ymax - ymin;
    bs = max(dx,dy);
    cx = 0.5*(xmin + xmax);
    cy = 0.5*(ymin + ymax);
    clmparams.xylim = [cx - 0.6*bs, cx + 0.6*bs, cy-0.6*bs, cy + 0.6*bs];
    
    clmparams.ndomain = geom_class.ndomain;
    if(isfield(geom_class,'rn'))
       clmparams.rn = geom_class.rn;
    else
       clmparams.rn = ones(geom_class.ndomain,1);
    end
    
    if(isfield(geom_class,'lambda'))
       clmparams.lambda = geom_class.lambda;
    else
        clmparams.lambda = 0.38;
    end
    ndomain = geom_class.ndomain;
    
    rn = clmparams.rn;
    lambda = clmparams.lambda;
    
    if(isfield(geom_class,'mode'))
        clmparams.mode = geom_class.mode;
    else
        clmparams.mode = 'te';
    end
    
    [~,nverts] = size(geom_class.verts);
    clmparams.verts = zeros(2,nverts+2);
    clmparams.verts(:,1:nverts) = geom_class.verts;
    clmparams.verts(1,nverts+1) = clmparams.xylim(1);
    clmparams.verts(1,nverts+2) = clmparams.xylim(2);
    verts = clmparams.verts;
    
    if(strcmpi(clmparams.mode,'te'))
        coef = ones(clmparams.ndomain,1);
    elseif(strcmpi(clmparams.mode,'tm'))
        coef = rn.^2;
    else
        fprintf('Invalid mode specification, reverting to te')
        coef = ones(geom_class.ndomain,1);
    end
    clmparams.coef = coef;
    
    % Determine the extent to go out on either side
    
    clmparams.lvert = nverts+1;
    clmparams.rvert = nverts+2;
    
    
    xmin = clmparams.xylim(1);
    xmax = clmparams.xylim(2);
    
    is_inf = zeros(1,ndomain);
    k = rn/lambda*2*pi;
    kinf = inf*ones(1,ndomain);
    for i=1:ndomain
        if(geom_class.regions{i}.is_inf ~= 0)
            is_inf(i) = geom_class.regions{i}.is_inf;
            kinf(i) = k(i); 
        end
    end
    
    
    if(isfield(geom_class,'ngr'))
        clmparams.ngr = geom_class.ngr;
    else
        clmparams.ngr = 300;
    end
    bs = xmax-xmin;
    krmax = max(real(k(:)));
    ngruse = ceil(bs*krmax/pi);
    
    clmparams.ngr = max(clmparams.ngr,ngruse);
    
    c1 = log(1.0/eps)/min(kinf)*2*pi;
    c2 = c1/2;
    
    L(1) = 6*c2-xmin;
    L(2) = 6*c2+xmax;
    
    % Determine the number of curves actually used
    % There are two additional curves corresponding to
    % left and right intervals and moreover additional 
    % curves get created if there is a polygon with
    % more than one edge each edge needs to be a separate
    % curve
    
    ncurve = geom_class.ncurve;
    ncurve_use = 4;
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
    cparams{1}.tb = xmin;
    cparams{1}.curvetype = 0;
    cparams{1}.ilr = 1;
    cparams{1}.ifclosed = false;
    cpars{1}.L = L(1); cpars{1}.c1 = c1; cpars{1}.c2 = c2;
    cpars{1}.lvert = clmparams.lvert;
    
    cparams{2}.ta = xmax;
    cparams{2}.tb = L(2);
    cparams{2}.curvetype = 0;
    cparams{2}.ilr = 2;
    cparams{2}.ifclosed = false;
    cpars{2}.L = L(2); cpars{2}.c1 = c1; cpars{2}.c2 = c2;
    cpars{2}.rvert = clmparams.rvert;
    
    cparams{3}.ta = 0;
    cparams{3}.tb = 1;
    cparams{3}.curvetype = 1;
    cparams{3}.ifclosed = false;
    cpars{3}.v0 = clmparams.verts(:,nverts+1);
    cpars{3}.v0ind = nverts+1;
    cpars{3}.v1 = clmparams.verts(:,geom_class.lvert);
    cpars{3}.v1ind = geom_class.lvert;
    
    cparams{4}.ta = 0;
    cparams{4}.tb = 1;
    cparams{4}.curvetype = 1;
    cparams{4}.ifclosed = false;
    cpars{4}.v0 = clmparams.verts(:,geom_class.rvert);
    cpars{4}.v0ind = geom_class.rvert;
    cpars{4}.v1 = clmparams.verts(:,nverts+2);
    cpars{4}.v1ind = nverts+2;
    
    
    icurve_cur = 5;
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
            if(isfield(geom_class.curves{i},'theta'))
               thet = geom_class.curves{i}.theta;
            else
               thet = pi/2;
            end
            cparams{icurve_cur}.ta = -thet/2;
            cparams{icurve_cur}.tb = thet/2;
            cparams{icurve_cur}.ifclosed = false;
            cparams{icurve_cur}.curvetype = ctype;
            cpars{icurve_cur}.v0 = verts(:,geom_class.curves{i}.vert_list(1));
            cpars{icurve_cur}.v0ind = geom_class.curves{i}.vert_list(1);
            cpars{icurve_cur}.v1 = verts(:,geom_class.curves{i}.vert_list(2));
            cpars{icurve_cur}.v1ind = geom_class.curves{i}.vert_list(2);
            cpars{icurve_cur}.theta = thet;
            if(isfield(geom_class.curves{i},'ifconvex'))
               cpars{icurve_cur}.ifconvex = geom_class.curves{i}.ifconvex;
            else
               cpars{icurve_cur}.ifconvex = 0;
            end
            icurve_cur = icurve_cur+1;
        elseif(ctype == 3)
            nwig = 1;
            if(isfield(geom_class.curves{i},'nwiggles'))
              nwig = geom_class.curves{i}.nwiggles;
            end
            amp = 1.0;
            if(isfield(geom_class.curves{i},'amplitude'))
               amp = geom_class.curves{i}.amplitude;
            end
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
        for j=1:length(geom_class.regions{i}.icurve_list)
            if(geom_class.regions{i}.icurve_list(j) < 0)
                clist{i} = [clist{i} -icurve_list{abs(geom_class.regions{i}.icurve_list(j))}]; 
            else
                clist{i} = [clist{i} icurve_list{abs(geom_class.regions{i}.icurve_list(j))}]; 
            end
        end 
        if(geom_class.regions{i}.is_inf == 1)
            clist{i} = [1 3 clist{i} 4 2];  
        elseif(geom_class.regions{i}.is_inf == -1)
            clist{i} = [-2 -4 clist{i} -3 -1];
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

    clmparams.issymmetric = 0;
    if(isfield(geom_class,'issymmetric'))
       clmparams.issymmetric = geom_class.issymmetric
    end
    clmparams.corners = corners;
    clmparams.ncorner = ncorner;
    
    % Determine nch
    nch = zeros(1,ncurve_use);
    fcurve = cell(1,ncurve_use);

    n0 = 12;
    fac = 1.2/2/pi;
    ppw = 5;
    if(isfield(geom_class,'ppw'))
        ppw = geom_class.ppw;
    end

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
      kmax= max(abs(k1(i)),abs(k2(i)));

      nch(i) = round(ppw*L*kmax/32/pi) + n0;
      %nch(i) = round(fac*L/lambda) + n0;
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
        if(isfield(geom_class.regions{i},'src_in'))
            src_in(:,i) = geom_class.regions{i}.src_in(:);
            issrc_in(i) = true;
        end
    end
    clmparams.src_in = src_in;
    clmparams.issrc_in = issrc_in;
    
    
    
    
end
