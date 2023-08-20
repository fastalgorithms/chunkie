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
                vmat = squeeze(allvs(:,:,ii,jj,kk));
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
end
