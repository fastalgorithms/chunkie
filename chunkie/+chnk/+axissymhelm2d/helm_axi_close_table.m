function [kerns,kernsdk,kernsda,kernsdkk,kernsdak,kernsdaa] ...
        = helm_axi_close_table(r0s,alphs,ifun,htables)
    int   = zeros(size(alphs(:)));
    kerns = int;
    kernsdk = int;
    kernsda = int;
    kernsdkk= int;
    kernsdak= int;
    kernsdaa= int;

    alphsuse = alphs(:);
    r0suse = r0s(:);
    
    iks = ceil(r0suse/pi);
    ias = ceil(-log(alphsuse*5)/log(2));
    

    iindin = find((iks<152).*(ias<111));
    
    r0sin = r0suse(iindin);
    alphsin = alphsuse(iindin);

    iksin = iks(iindin);
    iasin = ias(iindin);

    krel = (r0sin-pi*(iksin-1))*2/pi-1;
    arel = ((alphsin - 0.2*2.^(-iasin))./(0.2*2.^(-iasin))-0.5)*2;

    ncheb = htables.ncheb;
        
    tas = cos((0:(htables.ncheb-1)).'*acos(arel).');
    tks = cos((0:(htables.ncheb-1)).'*acos(krel).');

    naind = 111;
    nkind =  152;
    iindlin = (iksin - 1)*naind + iasin;
    vs = htables.allvs(:,:,iindlin,ifun);
    vdk = htables.alldk(:,:,iindlin,ifun);
    vda = htables.allda(:,:,iindlin,ifun);
    vdak = htables.alldak(:,:,iindlin,ifun);
    vdkk = htables.alldkk(:,:,iindlin,ifun);
    vdaa = htables.alldaa(:,:,iindlin,ifun);

    % 
    % 
    % vs = vs(:,:,iindlin);
    % vdk = vdk(:,:,iindlin);
    % vda = vda(:,:,iindlin);
    % vdak = vdak(:,:,iindlin);
    % vdaa = vdaa(:,:,iindlin);
    % vdkk = vdkk(:,:,iindlin);
    % 

    tkrep = repmat(tks,[1,1,ncheb]);
    tkrep = permute(tkrep,[3,1,2]);
    
    fvs = squeeze(dot(tkrep,vs,2));
    kern = dot(tas,fvs,1);
    kerns(iindin) = kern(:);
    

    fvdk = squeeze(dot(tkrep,vdk,2));
    kerndk = dot(tas,fvdk,1);
    kernsdk(iindin) = kerndk(:);

    fvda = squeeze(dot(tkrep,vda,2));
    kernda = dot(tas,fvda,1);
    kernsda(iindin) = kernda(:);

    fvdak = squeeze(dot(tkrep,vdak,2));
    kerndak = dot(tas,fvdak,1);
    kernsdak(iindin) = kerndak(:);

    fvdaa = squeeze(dot(tkrep,vdaa,2));
    kerndaa = dot(tas,fvdaa,1);
    kernsdaa(iindin) = kerndaa(:);


    fvdkk = squeeze(dot(tkrep,vdkk,2));
    kerndkk = dot(tas,fvdkk,1);
    kernsdkk(iindin) = kerndkk(:);

    kerns = reshape(kerns, size(alphs));
    kernsda = reshape(kernsda, size(alphs));
    kernsdk = reshape(kernsdk, size(alphs));
    kernsdkk = reshape(kernsdkk, size(alphs));
    kernsdak = reshape(kernsdak, size(alphs));
    kernsdaa = reshape(kernsdaa, size(alphs));

end

