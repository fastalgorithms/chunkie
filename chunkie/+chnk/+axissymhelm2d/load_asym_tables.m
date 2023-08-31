function asym_tables = load_asym_tables()
%CHNK.AXISSYMHELM2D.load_asym_tables loads the precomputed tables
% for axissymmetric Helmholtz kernels

    dirname = '+chnk/+axissymhelm2d/';
    vdats = {};
    nks = 152;
    for ii=1:nks
        load([dirname 'k2test',num2str(ii),'.mat']);
        v = Expression1;
        v = reshape(v,[12,12,3,111]);
        vdats{ii} = v;
    end

    allvs = zeros([size(vdats{1}),nks]);

    for ii=1:nks
        allvs(:,:,:,:,ii) = vdats{ii};
    end
    ncheb = 12;
    xcheb = (0:(ncheb-1))/(ncheb-1)*pi;

    [N,X]=meshgrid(0:(ncheb-1),xcheb);
    Tc2v = cos(N.*X);
    Tv2c = inv(Tc2v);
    Tfull = kron(Tv2c,Tv2c);
    Tc2vd = N.*sin(N.*X)./sqrt(1-cos(X).^2);
    Tc2vd(1,:) = (0:(ncheb-1)).^2;
    Tc2vd(end,:) = (-1).^(1:ncheb).*(0:(ncheb-1)).^2;
    Tc2cd = Tv2c*Tc2vd;

    xcheb = cos(xcheb)';


    allda = allvs;
    alldk = allvs;
    alldaa= allvs;
    alldak= allvs;
    alldkk= allvs;

    for kk=1:nks
        for jj=1:111
            for ii=1:3  
%                vmat = squeeze(allvs(:,:,ii,jj,kk));
                if (ii ==1)
                    vmat = squeeze(allvs(:,:,ii,jj,kk));
                else
                    vmat = squeeze(allvs(:,:,ii,jj,kk));
                    vmat = reshape(Tfull*vmat(:),[12,12]);
                    allvs(:,:,ii,jj,kk) = vmat;
                end
                ia = jj;
                vda = (2)^(ia-1)*20*Tc2cd*vmat;
                vdaa= (2)^(ia-1)*20*Tc2cd*vda;
                vdk = 2/(pi)*Tc2cd*vmat.';
                vdkk= (2)/pi*Tc2cd*vdk;
                vdak = 2/(pi)*Tc2cd*vda.';
                allda(:,:,ii,jj,kk) = vda;
                alldaa(:,:,ii,jj,kk)= vdaa;
                alldk(:,:,ii,jj,kk) = vdk.';
                alldak(:,:,ii,jj,kk)= vdak.';
                alldkk(:,:,ii,jj,kk)= vdkk.';
            end
        end
    end    

    naind = 111;
    nkind = 152;

    allvs = reshape(allvs, [ncheb, ncheb, 3, naind*nkind]);
    allda = reshape(allda, [ncheb, ncheb, 3, naind*nkind]);
    alldk = reshape(alldk, [ncheb, ncheb, 3, naind*nkind]);
    alldak = reshape(alldak, [ncheb, ncheb, 3, naind*nkind]);
    alldkk = reshape(alldkk, [ncheb, ncheb, 3, naind*nkind]);
    alldaa = reshape(alldaa, [ncheb, ncheb, 3, naind*nkind]);

    allvs = permute(allvs, [1,2,4,3]);
    allda = permute(allda, [1,2,4,3]);
    alldk = permute(alldk, [1,2,4,3]);
    alldak = permute(alldak, [1,2,4,3]);
    alldaa = permute(alldaa, [1,2,4,3]);
    alldkk = permute(alldkk, [1,2,4,3]);
    

    asym_tables =[];
    asym_tables.allvs = allvs;
    asym_tables.allda = allda;
    asym_tables.alldaa= alldaa;
    asym_tables.alldk = alldk;
    asym_tables.alldak= alldak;
    asym_tables.alldkk= alldkk;
    asym_tables.ncheb = ncheb;

    nlege = 500;
    [xlege,wlege,~,~] = lege.exps(nlege);
    xlege = (pi*(xlege+1)/2);
    wlege = wlege*pi/2;
    asym_tables.xlege = xlege;
    asym_tables.wlege = wlege;
    
   	nlege = 100;
    [xlege,wlege,~,~] = lege.exps(nlege);
    xlege = (pi*(xlege+1)/2);
    wlege = wlege*pi/2;
    asym_tables.xlege_mid = xlege;
    asym_tables.wlege_mid = wlege;
    
  	nlege = 50;
    [xlege,wlege,~,~] = lege.exps(nlege);
    xlege = (pi*(xlege+1)/2);
    wlege = wlege*pi/2;
    asym_tables.xlege_midnear = xlege;
    asym_tables.wlege_midnear = wlege;
    
  	nlege = 20;
    [xlege,wlege,~,~] = lege.exps(nlege);
    xlege = (pi*(xlege+1)/2);
    wlege = wlege*pi/2;
    asym_tables.xlege_midnearnear = xlege;
    asym_tables.wlege_midnearnear = wlege;
end
