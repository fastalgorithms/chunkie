
function geom_class = read_geom_clm8()
    geom_class = [];
    geom_class.xylim = [-15 15 -26 4];
    ndomain = 5;
    geom_class.ndomain = ndomain;
    ncurve = 11;
    geom_class.ncurve = ncurve;
    

    rn = ones(ndomain,1); 
  % rn(i) is the index of refraction of the ith domain
    rn(1) = 1.0;
    rn(2) = 1.34+1i*1e-2;
    rn(3) = 1.452;
    rn(4) = 1.348;
    rn(5) = 1.363;
    
    lambda = 0.38;
    geom_class.rn = rn;
    geom_class.lambda = lambda;
    geom_class.mode = 'te';
    
    a = 2.1426;
    
    verts = zeros(2,10);
    verts(1,1) = -a;
    verts(1,2) = a;
    
    verts(:,3) = [-0.72;-2.6061];
    verts(:,4) = [0.72;-2.6061];
    
    verts(:,7) = [-0.9;-2.8468];
    verts(:,8) = [0.9;-2.8468];
    
    verts(:,5) = [-0.09;-16.2665];
    verts(:,6) = [0.09;-16.2665];
    
    
    verts(1,9) = -5;
    verts(1,10) = 5;
    
    % two circular arcs for the center eye for now
    theta = zeros(1,3);
    % upper curve opening angle
    theta(1) = pi/1.3;
    theta(3) = pi/8;


    % curve #4 -- sine curve on [a,b]
    n4 = 3; A4 = 0.3;

    
    geom_class.verts = verts;

    curves = cell(1,ncurve);
    curves{1}.curvetype = 2;
    curves{2}.curvetype = 3;
    curves{3}.curvetype = 1;
    curves{4}.curvetype = 1;
    curves{5}.curvetype = 1;
    curves{6}.curvetype = 1;
    curves{7}.curvetype = 1;
    curves{8}.curvetype = 2;
    curves{9}.curvetype = 1;
    curves{10}.curvetype = 1;
    curves{11}.curvetype = 1;
    
    curves{1}.vert_list = [1 2];
    curves{2}.vert_list = [1 2];
    curves{3}.vert_list = [1 3];
    curves{4}.vert_list = [4 2];
    curves{5}.vert_list = [3 4];
    curves{6}.vert_list = [7 5];
    curves{7}.vert_list = [6 8];
    curves{8}.vert_list = [5 6];
    curves{9}.vert_list = [7 8];
    curves{10}.vert_list = [9,1];
    curves{11}.vert_list = [2,10];
    
    for i=1:ncurve
        curves{i}.curve_id = i;
    end
    curves{1}.theta = theta(1);
    curves{1}.ifconvex = 0;
    
    curves{8}.theta = theta(3);
    curves{8}.ifconvex = 2;
    
    curves{2}.nwiggles = n4;
    curves{2}.amplitude = A4;
    
    geom_class.curves = curves;
    
    
%   Define region info    
    clist = cell(1,ndomain);
    clist{1} = [10 -1 11];
    clist{2} = [-11 -4 -5 -3 -10 -7 -8 -6 9];
    clist{3} = [1,2];
    clist{4} = [-2 3 5 4];
    clist{5} = [-9 6 8 7];
    regions = cell(1,ndomain);
    src = zeros(2,ndomain);

    src(2,:) = [4.67,  -2.5e1, 0.7, -1.68, -4.68];

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
    geom_class.lvert = 9;
    geom_class.rvert = 10;
    

end