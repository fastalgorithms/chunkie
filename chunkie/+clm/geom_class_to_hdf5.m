function [] = geom_class_to_hdf5(geom_class,fname)
   if(exist(fname,'file')==2)
       fprintf('file already exists: overwriting\n');
   end
   fcpl = H5P.create('H5P_FILE_CREATE');
   fapl = H5P.create('H5P_FILE_ACCESS');
   fid = H5F.create(fname,'H5F_ACC_TRUNC',fcpl,fapl);
   
   h5writeatt(fname,'/','ndomain',geom_class.ndomain);
   h5writeatt(fname,'/','rnr',real(geom_class.rn));
   h5writeatt(fname,'/','rni',imag(geom_class.rn));
   h5writeatt(fname,'/','mode',geom_class.mode);
   h5writeatt(fname,'/','ncurve',geom_class.ncurve);
   h5writeatt(fname,'/','verts',geom_class.verts);
   h5writeatt(fname,'/','xylim',geom_class.xylim);
   h5writeatt(fname,'/','lambda',geom_class.lambda);
   h5writeatt(fname,'/','lvert',geom_class.lvert);
   h5writeatt(fname,'/','rvert',geom_class.rvert);
   plist = 'H5P_DEFAULT';
   gid = H5G.create(fid,'curves',plist,plist,plist);
   ncurve = geom_class.ncurve;
   for i=1:ncurve
       str1 = int2str(i);
       gid1 = H5G.create(gid,str1,plist,plist,plist);
       str0 = ['/curves/' str1];
       h5writeatt(fname,str0,'curvetype',geom_class.curves{i}.curvetype);
       h5writeatt(fname,str0,'vert_list',geom_class.curves{i}.vert_list);
       h5writeatt(fname,str0,'curve_id',geom_class.curves{i}.curve_id);
       ctype = geom_class.curves{i}.curvetype;
       if(ctype == 2)
           h5writeatt(fname,str0,'theta',geom_class.curves{i}.theta);
           h5writeatt(fname,str0,'ifconvex',geom_class.curves{i}.ifconvex);
       elseif(ctype == 3)
           h5writeatt(fname,str0,'nwiggles',geom_class.curves{i}.nwiggles);
           h5writeatt(fname,str0,'amplitude',geom_class.curves{i}.amplitude);
       end
       H5G.close(gid1);
   end
   H5G.close(gid);
   gid = H5G.create(fid,'regions',plist,plist,plist);
   ndomain = geom_class.ndomain;
   for i=1:ndomain
       str1 = int2str(i);
       gid1 = H5G.create(gid,str1,plist,plist,plist);
       str0 = ['/regions/' str1];
       h5writeatt(fname,str0,'region_id',geom_class.regions{i}.region_id);
       h5writeatt(fname,str0,'icurve_list',geom_class.regions{i}.icurve_list);
       h5writeatt(fname,str0,'is_inf',geom_class.regions{i}.is_inf);
       h5writeatt(fname,str0,'src_in',geom_class.regions{i}.src_in);
       H5G.close(gid1);
   end
   H5F.close(fid);
end