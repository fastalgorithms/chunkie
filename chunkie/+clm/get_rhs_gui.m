function rhs = get_rhs_gui(chnkr,clmparams,np,alpha1,alpha2,opts)

    ngl = 16;
    if isfield(clmparams,'k')
        k = clmparams.k;
    end

    if isfield(clmparams,'c')
      c = clmparams.c;
    end
    if isfield(clmparams,'k1')
      k1 = clmparams.k1;
    end
    if isfield(clmparams,'k2')
      k2 = clmparams.k2;
    end
    if isfield(clmparams,'coef')
      coef = clmparams.coef;
    end
    
    if isfield(clmparams,'ncurve')
      ncurve = clmparams.ncurve;
    end
    if isfield(clmparams,'ndomain')
      ndomain = clmparams.ndomain;
    end
    

    if isfield(clmparams, 'nch')
      nch = clmparams.nch;
    end

    if isfield(clmparams, 'src')
      src = clmparams.src;
    end

    rhs = zeros(2*np,1);
    if(opts.itype == 1)
% construct artificial boundary data for testing purpose
    

        for i=1:ncurve
          d1 = c(1,i); % interior domain index
          d2 = c(2,i); % exterior domain index

          j1 = d1 + 1; % src index for the interior domain
          if j1 > ndomain
            j1 = j1 - ndomain;
          end
          j2 = d2 + 1; % src index for the exterior domain
          if j2 > ndomain
            j2 = j2 - ndomain;
          end

          c1 = coef(d1);
          c2 = coef(d2);

          ind1 = sum(nch(1:i-1))*ngl*2+(1:2:2*nch(i)*ngl);
          ind2 = sum(nch(1:i-1))*ngl*2+(2:2:2*nch(i)*ngl);

          targnorm = chnkr(i).n;
          nx = targnorm(1,:); nx=nx.';
          ny = targnorm(2,:); ny=ny.';

          [u1,gradu1]=chnk.helm2d.green(k1(i),src(:,j1),chnkr(i).r(:,:));
          du1dn = gradu1(:,1).*nx + gradu1(:,2).*ny;

          [u2,gradu2]=chnk.helm2d.green(k2(i),src(:,j2),chnkr(i).r(:,:));
          du2dn = gradu2(:,1).*nx + gradu2(:,2).*ny;

          rhs(ind1) = -alpha1(i)*(u1-u2);
          rhs(ind2) =  alpha2(i)*(1/c1*du1dn-1/c2*du2dn);
        end

    elseif (opts.itype == 2)
        alpha = opts.alpha;
        for i=1:ncurve
            d1 = c(1,i); % interior domain index
            d2 = c(2,i); % exterior domain index

            c1 = coef(d1);
            c2 = coef(d2);

            ind1 = sum(nch(1:i-1))*ngl*2+(1:2:2*nch(i)*ngl);
            ind2 = sum(nch(1:i-1))*ngl*2+(2:2:2*nch(i)*ngl);

            targnorm = chnkr(i).n;
            nx = targnorm(1,:,:); nx=nx(:);
            ny = targnorm(2,:,:); ny=ny(:);

            if d1==1 || d1 == 2
                [u1,gradu1]=clm.planewavetotal(k(1),alpha,k(2),chnkr(i).r,d1,coef);
                du1dn = gradu1(:,1).*nx + gradu1(:,2).*ny;
                rhs(ind1) = rhs(ind1) + alpha1(i)*u1(:);
                rhs(ind2) = rhs(ind2) - alpha2(i)/c1*du1dn(:);
            end

            if d2==1 || d2==2
                [u2,gradu2]=clm.planewavetotal(k(1),alpha,k(2),chnkr(i).r,d2,coef);
                du2dn = gradu2(:,1).*nx + gradu2(:,2).*ny;
                rhs(ind1) = rhs(ind1) - alpha1(i)*u2(:);
                rhs(ind2) = rhs(ind2) + alpha2(i)/c2*du2dn(:);
            end
        end
    end
    rhs = rhs(:);
end
