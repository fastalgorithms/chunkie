function [geom_class,ier] = hdf5_to_geom_class(fname)
   geom_class = [];
   ier = 0;
   try
      geom_class.ndomain = h5readatt(fname,'/','ndomain');
   catch
      fprintf('ndomain attribute missing.\n This error is fatal.\n Exiting routine.\n\n\n');
      ier = 1;
      return
   end
   ndomain = geom_class.ndomain;
   try
      geom_class.ncurve = h5readatt(fname,'/','ncurve');
   catch
      fprintf('ncurve attribute missing.\n This error is fatal.\n Exiting routine.\n\n\n');
      ier = 2;
      return
   end
   ncurve = geom_class.ncurve;
   try
       geom_class.verts = h5readatt(fname,'/','verts');
   catch
       fprintf('verts attribute missing. \n This error is fatal.\n Beed at least 2 vertices.\n Exiting routine.\n\n\n');
       ier = 3;
       return
   end
   
   try
      rnr = h5readatt(fname,'/','rnr');
   catch
      fprintf('real part of refractive index not set. Setting to default value of 1.\n\n\n');
      rnr = ones(ndomain,1);
   end
   try
      rni = h5readatt(fname,'/','rni');
   catch
      fprintf('imaginary part of refractive index not set. Setting to default value of 0.\n\n\n');
      rni = zeros(ndomain,1);
   end
   rn = rnr + 1j*rni;
   geom_class.rn = rn;
   try
      geom_class.mode = h5readatt(fname,'/','mode');
   catch
      fprintf('mode not set. Setting to "te" mode.\n\n\n');
      geom_class.mode = 'te'; 
   end
   try
      geom_class.xylim = h5readatt(fname,'/','xylim');
   catch
       fprintf('xylim not set.\n\n\n'); 
       xmin = min(geom_class.verts(1,:));
       xmax = max(geom_class.verts(1,:));
       ymin = min(geom_class.verts(2,:));
       ymax = max(geom_class.verts(2,:));
       dx = xmax - xmin;
       dy = ymax - ymin;
       bs = max(dx,dy);
       cx = (xmin + xmax)*0.5;
       cy = (ymin + ymax)*0.5;
       geom_class.xylim = [cx - bs, cx + bs, cy - bs, cy+bs];
      
   end
   
   try
      geom_class.ngr = h5readatt(fname,'/','ngr');
   catch
       fprintf('ngr not set. Setting to default value of 300\n\n\n'); 
       geom_class.ngr = 300;
   end
   
   try
      geom_class.lambda = h5readatt(fname,'/','lambda');
   catch
       fprintf('lambda not set.\n\n\n');
       geom_class.lambda = 0.38;
       
   end
   
   
   try
      geom_class.lvert = h5readatt(fname,'/','lvert');
   catch
       fprintf('lvert not set.\n This error is fatal.\n exiting.\n\n\n');
       ier = 4;
       return
   end
   
   try
      geom_class.rvert = h5readatt(fname,'/','rvert');
   catch
       fprintf('rvert not set.\n This error is fatal.\n exiting.\n\n\n');
       ier = 5;
       return
   end
   curves = cell(1,ncurve);
   for i=1:ncurve
       str1 = int2str(i);
       str0 = ['/curves/' str1];
       try
          curves{i}.curvetype = h5readatt(fname,str0,'curvetype');
       catch
          fprintf('curve type not specified.\n This error is fatal.\n exiting.\n\n\n');
          ier = 6;
          return;
       end
       try
          curves{i}.curve_id = h5readatt(fname,str0,'curve_id');
       catch
          fprintf('curve id not specified.\n This error is fatal.\n exiting.\n\n\n');
          ier = 7;
          return;
       end
       ctype = curves{i}.curvetype;
       if (ctype == 1 || ctype == 2 || ctype == 3)
           try
              curves{i}.vert_list = h5readatt(fname,str0,'vert_list')';
           catch
               fprintf('vert list not specified for curve type 1,2, or 3.\n This error is fatal. Exiting.\n\n\n');
               ier = 8;
               return;
           end
       end
       
       if(ctype == 2)
           try
               curves{i}.theta = h5readatt(fname,str0,'theta');
           catch
               fprintf('missing value of theta for curve type 2.\n Setting it to default value pi/4.\n\n\n');
               curves{i}.theta = pi/2;
           end
           try
               curves{i}.ifconvex = h5readatt(fname,str0,'ifconvex');
           catch
               fprintf('missing value of ifconvex for curve type 2.\n Setting it to default value 0.\n\n\n');
               curves{i}.ifconvex = 0;
           end
       elseif(ctype == 3)
           try
               curves{i}.nwiggles = h5readatt(fname,str0,'nwiggles');
           catch
               fprintf('missing value of nwiggles for curve type 3.\n Setting it to default value 2.\n\n\n');
               curves{i}.nwiggles = 2;
           end
           try
               curves{i}.amplitude = h5readatt(fname,str0,'amplitude');
           catch
               fprintf('missing value of amplitude for curve type 3.\n Setting it to default value 1.\n\n\n');
               curves{i}.amplitude = 1.0;
           end
       end
   end
   geom_class.curves = curves;
   regions = cell(1,ndomain);
   
   for i=1:ndomain
       str1 = int2str(i);
       str0 = ['/regions/' str1];
       try
           regions{i}.region_id = h5readatt(fname,str0,'region_id');
       catch
          fprintf('region id not specified.\n This error is fatal.\n exiting.\n\n\n');
          ier = 9;
          return;
       end
       
       try
           regions{i}.icurve_list = h5readatt(fname,str0,'icurve_list')';
       catch
          fprintf('icurve_list not specified.\n This error is fatal.\n exiting.\n\n\n');
          ier = 10;
          return;
       end
       
       try
           regions{i}.is_inf = h5readatt(fname,str0,'is_inf');
       catch
           fprintf('is_inf not specified.\n Setting it to default value of 0');
           regions{i}.is_inf = 0;
       end
       try
           regions{i}.src_in = h5readatt(fname,str0,'src_in');
       catch
       end
   end
   geom_class.regions = regions;
end