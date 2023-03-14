function [kerns,varargout] = transmission_helper(chnkobj,ks,cs,coefs,varargin)
%CHNK.HELM2D.TRANSMISSION_HELPER builds the matrix of kernels required to call
% chunkermat and generate boundary data which is appropriately scaled for the linear
% system if requested.
%
% Syntax: kerns = transmission_helper(chnkrs,ks,cs,coefs), returns a matrix of kernels
%         [kerns,bdrydata] = transmission_helper(chnkrs,ks,cs,coefs,opts), returns
%            a matrix of kernels and the boundary data due to point sources or
%            waves as prescribed by opts
%
% Note: Here we always assume that the transmission problem is either in
%    the 'TE' or 'TM' polarization, which implies that the boundary conditions
%    are of the form
%    [u] = f, [coef du/dn] = g
%
%    where [u] denotes the jump in u across an interface. 
%    To obtain the TE polarization, coef \equiv 1, and for the TM polarization
%    coef should be set to the square of the refractive index of the medium.
%
%    These set of transmission problems do not cover all possibilities, and
%    in particular exclude transmission problems with boundary conditions
%    of the form
%    [\beta u] = f. 
%
% Input:
%   chnkrs - cell array of chunker objects describing the geometry (1,ncurve)
%   ks - wavenumbers in different regions of teh plane
%   cs - cs(1,i), cs(2,i) denote the region number in the direction of the normal
%          and opposite of the direction of the normal for curve chnkrs(i)
%   coefs - denotes the scaling parameter for the jump in the neumann data, depends
%            only on the region
%   opts - optional structure for defining the boundary condition
%      opts.bdry_data_type - 'pw' or 'point sources'
%      opts.bdry_data_type = pw, then opts.dir should be set to incident
%      angle, and opts.exposed_curves should be set to the curves for which
%          incident field is to be set
%      if opts.bdry_data_type = 'point sources' then opts.sources should be a cell
%      array of size (1,ncurve), with  opts.sources{i}.locs (2,m) locations 
%      of charges which generate the data for region i, and .charges are the
%      corresponding charge strengths
%           
%
    
    
    if nargin == 4
      opts = [];
        
    elseif nargin == 5
      opts = varargin{1};
    else
       fprintf('Invalid number of input arguments\n');
       fprintf('Returning without computing anything\n');
       % this behavior might need to be fixed
       kerns = [];
       varargout{2:nargout} = [];
       return;  
    end
    
    
    if (class(chnkobj) == "chunker")
        chnkrs = chnkobj;
    elseif(class(chnkobj) == "chunkgraph")
        chnkrs = chnkobj.echnks;
    else
        msg = "Unsupported object in chunkermat";
        error(msg)
    end

    
    bdry_data_type = 'pw';
    direction = 0;
    ncurve = length(chnkrs);
    sources = cell(1,ncurve);
    charges = cell(1,ncurve);
    exposed_curves = ones(1,ncurve);
    for i=1:ncurve
        sources{i} = zeros(2,1);
        charges{i} = 0 + 1j*0;
    end
    
    if(isfield(opts,'bdry_data_type')) 
        bdry_data_type = opts.bdry_data_type;
    end
    if(isfield(opts,'dir'))
        direction = opts.dir;
    end
    if(isfield(opts,'sources'))
        sources = opts.sources;
    end
    if(isfield(opts,'charges'))
        charges = opts.charges;
    end
    if(isfield(opts,'exposed_curves'))
        exposed_curves = opts.exposed_curves;
    end
     
    
     
    
    kerns = cell(ncurve,ncurve);
        
    k1 = ks(cs(1,:));
    k2 = ks(cs(2,:));
    
    alpha1 = 2./(coefs(cs(1,:)) + coefs(cs(2,:)));
    alpha2 = 2./(1./coefs(cs(1,:)) + 1./coefs(cs(2,:)));
    
    
% First build system matrix without any corner corrections
    cc1 = zeros(2,2);
    cc2 = zeros(2,2);  
    for i=1:ncurve
        
        c1 = coefs(cs(1,i));
        c2 = coefs(cs(2,i));
        cc1(1,1) = -alpha1(i)*c1;
        cc2(1,1) = -alpha1(i)*c2;
        
        cc1(1,2) = -alpha1(i);
        cc2(1,2) = -alpha1(i);
        
        cc1(2,1) = alpha2(i);
        cc2(2,1) = alpha2(i);
        
        cc1(2,2) = alpha2(i)/c1;
        cc2(2,2) = alpha2(i)/c2;
        
        for j=1:ncurve
            kerns{i,j} = @(s,t) -(chnk.helm2d.kern(k1(i),s,t,'all',cc1)- ...
                 chnk.helm2d.kern(k2(i),s,t,'all',cc2));
        end  
    end
    
    
    
    nchs = zeros(1,ncurve);
    for i=1:ncurve
        nchs(i) = chnkrs(i).nch; 
    end
    
    ngl = chnkrs(1).k;
    np = sum(nchs(1:ncurve))*ngl;
    
    
    if nargout > 1
        bdry_data = complex(zeros(2*np,1));
        if(strcmpi(bdry_data_type,'pw'))
            for i=1:ncurve
                d1 = cs(1,i);
                d2 = cs(2,i);
                c1 = coefs(d1);
                c2 = coefs(d2);
                ind1 = sum(nchs(1:i-1))*ngl*2+(1:2:2*nchs(i)*ngl);
                ind2 = sum(nchs(1:i-1))*ngl*2+(1:2:2*nchs(i)*ngl);
                targnorm = chnkrs(i).n;
                nx = targnorm(1,:); nx = nx(:);
                ny = targnorm(2,:); ny = ny(:);
                if(exposed_curves(i))
                    
                    if(exposed_curves(i) < 0) 
                        k = k2(i);
                        a1 = alpha1(i);
                        a2 = -alpha2(i)/c2;
                    else
                        k = k1(i);
                        a1 = -alpha1(i);
                        a2 = alpha2(i)/c1;
                        
                    end
                    x = chnkrs(i).r(1,:); 
                    x = x(:);
                    y = chnkrs(i).r(2,:);
                    y = y(:);
                    ct = cos(direction);
                    st = sin(direction);
                    u = exp(1i*k*(x*ct + y*st));
                    dudn = 1i*k*(ct*nx+st*ny).*u;
                    bdry_data(ind1) = a1*u;
                    bdry_data(ind2) = a2*dudn;
                    
                end
            end
            
        elseif(strcmpi(bdry_data_type,'point sources'))
            
            for i=1:ncurve
                d1 = cs(1,i);
                d2 = cs(2,i);
                c1 = coefs(d1);
                c2 = coefs(d2);
                
                src1 = sources{d1};
                charges1 = charges{d1};
                charges1 = charges1(:);
                src2 = sources{d2};
                charges2 = charges{d2};
                charges2 = charges2(:);
                
                ind1 = sum(nchs(1:i-1))*ngl*2+(1:2:2*nchs(i)*ngl);
                ind2 = sum(nchs(1:i-1))*ngl*2+(2:2:2*nchs(i)*ngl);
                
                targnorm = chnkrs(i).n;
                nx = targnorm(1,:); nx = nx.';
                ny = targnorm(2,:); ny = ny.';
                [u1,gradu1] = chnk.helm2d.green(k1(i),src1,chnkrs(i).r(:,:));
                u1 = u1*charges1;
                gradu1(:,:,1) = gradu1(:,:,1)*charges1;
                gradu1(:,:,2) = gradu1(:,:,2)*charges1;
                dudn1 = squeeze(gradu1(:,1)).*nx + squeeze(gradu1(:,2)).*ny;
                
                [u2,gradu2] = chnk.helm2d.green(k2(i),src2,chnkrs(i).r(:,:));
                u2 = u2*charges2;
                gradu2(:,:,1) = gradu2(:,:,1)*charges2;
                gradu2(:,:,2) = gradu2(:,:,2)*charges2;
                dudn2 = squeeze(gradu2(:,1)).*nx + squeeze(gradu2(:,2)).*ny;
                
                bdry_data(ind1) = alpha1(i)*(u1-u2);
                bdry_data(ind2) = -alpha2(i)*(1/c1*dudn1 - 1/c2*dudn2);
            end
            
        end
        varargout{1} = bdry_data;
        
    end
    

end
