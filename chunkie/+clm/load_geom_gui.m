function geom_class = load_geom_gui(fname)
    [~,~,ext] = fileparts(fname);
    if(strcmpi(ext,'.h5'))
        disp(fname)
        [geom_class,ier] = clm.hdf5_to_geom_class(fname);
        if(ier ~= 0)
            fprintf('something terrible happened. Try different geometry file\n\n');
        end
    elseif(strcmpi(ext,'.mat'))
        a = load(fname);
        geom_class = a.geom_class;
    end
end