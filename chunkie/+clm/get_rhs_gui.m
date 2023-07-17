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

    if isfield(clmparams, 'src_in')
      src = clmparams.src_in;
    end

    
    if isfield(clmparams, 'src')
      src = clmparams.src_in;
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
        idomup = find(clmparams.is_inf == 1);
        idomdown = find(clmparams.is_inf == -1);
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

            if d1==idomup || d1 == idomdown
                y = chnkr(i).r(2,:);
                if(norm(y)>1e-16)
                    [u1,gradu1]=clm.planewavetotal_gui(k(idomup),alpha,k(idomdown),chnkr(i).r,clmparams.is_inf(d1),idomup,idomdown,coef);
                    du1dn = gradu1(:,1).*nx + gradu1(:,2).*ny;
                    rhs(ind1) = rhs(ind1) + alpha1(i)*u1(:);
                    rhs(ind2) = rhs(ind2) - alpha2(i)/c1*du1dn(:);
                
                    fprintf('maxu1=%d   icurve=%d\n',max(abs(u1)),i);
                end
            end

            if d2==idomup || d2==idomdown
                y = chnkr(i).r(2,:);
                if(norm(y)>1e-16)
                    [u2,gradu2]=clm.planewavetotal_gui(k(idomup),alpha,k(idomdown),chnkr(i).r,clmparams.is_inf(d2),idomup,idomdown,coef);
                    du2dn = gradu2(:,1).*nx + gradu2(:,2).*ny;
                    rhs(ind1) = rhs(ind1) - alpha1(i)*u2(:);
                    rhs(ind2) = rhs(ind2) + alpha2(i)/c2*du2dn(:);
                    fprintf('maxu2=%d   icurve=%d\n',max(abs(u2)),i);
                end
                
            end
        end
    end
    rhs = rhs(:);
end
