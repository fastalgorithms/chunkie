function asym_tables = load_asym_tables()
%CHNK.AXISSYMHELM2D.load_asym_tables loads the precomputed tables
% for axissymmetric Helmholtz kernels
    p = mfilename('fullpath');
    dirname = dir([p ,'.m']).folder;
    fname = [dirname '/asym_helm_data.mat'];
    load(fname,'allvs');
    msizes = size(allvs);
    ncheb = msizes(1);
    nks = msizes(5);
    nas = msizes(4);
    nkers = msizes(3);
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
    alldaa = allvs;
    alldak = allvs;
    alldkk = allvs;

    for kk=1:nks
        for jj=1:nas
            for ii=1:nkers
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
