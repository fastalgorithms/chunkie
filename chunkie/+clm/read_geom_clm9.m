
function geom_class = read_geom_clm9()
    geom_class = [];
    geom_class.xylim = [-4 4 -4 4];
    ndomain = 4;
    geom_class.ndomain = ndomain;
    ncurve = 5;
    geom_class.ncurve = ncurve;
    

    rn = ones(ndomain,1); 
  % rn(i) is the index of refraction of the ith domain
    rn(1) = 1.0;
    rn(2) = 1.34+1i*1e-2;
    rn(3) = 1.452;
    rn(4) = 1.348;
    
    
    lambda = 0.38;
    geom_class.rn = rn;
    geom_class.lambda = lambda;
    geom_class.mode = 'te';
    
    a = 2.1426;
    
    verts = zeros(2,4);
    verts(1,1) = -a;
    verts(1,2) = a;
    
    verts(:,3) = [-0.72;-2.6061];
    verts(:,4) = [0.72;-2.6061];
    
    
    % two circular arcs for the center eye for now
    theta = zeros(1,1);
    % upper curve opening angle
    theta(1) = pi/1.3;
    


    % curve #4 -- sine curve on [a,b]
    n4 = 3; A4 = 0.3;

    
    geom_class.verts = verts;

    curves = cell(1,ncurve);
    curves{1}.curvetype = 2;
    curves{2}.curvetype = 3;
    curves{3}.curvetype = 1;
    curves{4}.curvetype = 1;
    curves{5}.curvetype = 1;
    
    curves{1}.vert_list = [1 2];
    curves{2}.vert_list = [1 2];
    curves{3}.vert_list = [1 3];
    curves{4}.vert_list = [4 2];
    curves{5}.vert_list = [3 4];
    
    
    for i=1:ncurve
        curves{i}.curve_id = i;
    end
    curves{1}.theta = theta(1);
    curves{1}.ifconvex = 0;
        
    curves{2}.nwiggles = n4;
    curves{2}.amplitude = A4;
    
    geom_class.curves = curves;
    
    
%   Define region info    
    clist = cell(1,ndomain);
    clist{1} = [-1];
    clist{2} = [-4 -5 -3];
    clist{3} = [1,2];
    clist{4} = [-2 3 5 4];
    
    regions = cell(1,ndomain);
    src = zeros(2,ndomain);

    src(2,:) = [4.67,  -2.5e1, 0.7, -1.68];

    for i=1:ndomain
        regions{i}.region_id = 1;
        regions{i}.icurve_list = clist{i};
        regions{i}.is_inf = 0;
        if(i == 1)
            regions{i}.is_inf = 1;
        elseif (i==2)
            regions{i}.is_inf = -1;
        end
        regions{i}.src_in = src(:,i);
    end
    geom_class.regions = regions;
    geom_class.lvert = 1;
    geom_class.rvert = 2;
    

end