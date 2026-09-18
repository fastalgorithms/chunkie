function [qm, qmd, qmdd] = qleg_half_miller_vec(t, m)

    %
    % Uses Miller's Algorithm to compute the half Legendre functions of the
    % second kind
    %
    % Inputs
    % t: chi - 1
    % m: all the modes we want to compute
    %
    % Outputs
    % qm: [Q_{-1/2},...,Q_{m-1/2}]
    % qmd: [dQ_{-1/2},...,dQ_{m-1/2}]
    % qmdd: [d^2 Q_{-1/2},...,d^2 Q_{m-1/2}]
    %
    
    % compute chi
    chi = t+1;
    
    % check if we can run the forward recurrence
    if (m > 12307)
        mask1 = 1.00000005d0 <= (chi < 1.005);
    elseif (m > 4380)
        mask_forward = 1.0000005d0 <= (chi < 1.005);
    elseif (m > 1438)
        mask_forward = 1.000005d0 <= (chi < 1.005);
    elseif (m > 503)
        mask_forward = 1.00005d0 <= (chi < 1.005);
    elseif (m > 163)
        mask_forward = 1.0005d0 <= (chi < 1.005);
    else
        mask_forward = chi < 1.005; 
    end
    mask_back = ~mask_forward; % backward reccurence mask

    qm = zeros(m+1,size(t,2),size(t,3));
    qmd = zeros(m+1,size(t,2),size(t,3));
    qmdd = zeros(m+1,size(t,2),size(t,3));
    
    % run the forward reccurence
    if (sum(mask_forward,'all') ~= 0)
        [qm(:,mask_forward),qmd(:,mask_forward),qmdd(:,mask_forward)] = forward_reccurence(t(mask_forward),m);
    end

    % run the backward reccurence
    if (sum(mask_back,'all') ~= 0)
        [qm(:,mask_back),qmd(:,mask_back),qmdd(:,mask_back)] = backward_reccurence(t(mask_back),m);
    end

end

function [qm,qmd,qmdd] = forward_reccurence(t, m)
    t = t(:);
    chi = t+1;
    
    % initialize the seeds for the Legendre functions
    [q0,q1,q0d] = chnk.axissymlap2d.qleg_half(t);
    q1d = (-q0 + chi.*q1)./(2*(chi+1).*t);
    
    qm = zeros(m+1,size(t,1));
    qmd = zeros(m+1,size(t,1));
    qmdd = zeros(m+1,size(t,1));

    qm(1,:) = q0;
    qm(2,:) = q1;
    qmd(1,:) = q0d';
    qmd(2,:) = q1d';

    for i = 1:(m-1)
        j = i+1;
        qm(2+i,:) = 4*(j-1)/(2*j-1)*chi'.*qm(1+i,:) - (2*j-3)/(2*j-1)*qm(i,:);
        qmd(2+i,:) = ((j-1/2)*qm(1+i,:) - (j-1/2)*chi'.*qm(2+i,:))./(-t'.^2-2*t');
    end

    for i = 1:m+1
        j = i-1;
        qmdd(i,:) = (-(j-1/2)*(j+1/2)*qm(i,:) + 2*chi'.*qmd(i,:))./(-t'.^2-2*t');
    end
end

function [qm,qmd,qmdd] = backward_reccurence(t, m)
    t = t(:);
    chi = t+1;

    % initialize the seeds for the Legendre functions
    [q0,~,q0d] = chnk.axissymlap2d.qleg_half(t);

    % Millers Algorithm:
    qm = zeros(m+1,size(t,1));
    qmd = zeros(m+1,size(t,1));
    qmdd = zeros(m+1,size(t,1));

    % 1. run the forward reccurence until it has blown up
    % recurrence intialization
    fprev = ones(size(t));
    fprevprev = zeros(size(t));
    
    % NOTE:
    % should maybe wait until second derivatives also blow up

    maxiter = 1000; % max number of extra terms
    upbound = 1.0e17; % blow up tolerance
    nterms = zeros(size(t)); % number of extra terms
    
    % list of recurrences that have blown up
    already_blowup = zeros(size(t)); 
    for i = m:maxiter
        j = i+1;
        f = 4*(j-1)/(2*j-1)*chi.*fprev - (2*j-3)/(2*j-1)*fprevprev;
        d = ((j-1/2)*fprevprev - (j-1/2)*chi.*fprev)./(-t.^2-2*t);

        % check which values blew up
        f_blowup = (abs(f) >= upbound);

        % check which derivatives blew up
        d_blowup = (abs(d) >= upbound);

        % check where both blew up
        both_blowup = f_blowup.*d_blowup;

        % check only the new ones that blew up
        new_blowup = (both_blowup > already_blowup);
        
        % set the nterms for the new ones that blew up
        nterms = nterms + (i+1)*new_blowup;
        
        % update the list of ones that have blown up
        already_blowup = already_blowup + new_blowup;

        fprevprev = fprev;
        fprev = f;
    end

    % replace zeros in nterms with maxiter
    maxedout = nterms == 0;
    nterms = nterms + maxiter*maxedout;

    % 2. run the backward reccurence
    % recurrence intialization
    fprev = ones(size(t));
    fprevprev = zeros(size(t));
    
    % run the backward reccurence from nterms to m
    nterm_max = max(nterms, [], 'all');
    for i = 1:(nterm_max-m+1)
        j = nterm_max-i+1+1;

        f = 4*(j-1)/(2*j-3)*chi.*fprev - (2*j-1)/(2*j-3)*fprevprev;

        % check if we should start the backwards reccurence yet
        should_compute = (j <= nterms);
        
        % only update if we have started the reccurence for that element
        fprevprev = fprev.*should_compute;
        fprev = f.*should_compute + fprev.*(~should_compute);
    end

    qm(m,:) = fprev';
    qm(m+1,:) = fprevprev';

    % run the backward reccurence from m to 1
    for i = 1:(m-1)
        j = m-1-i+1+1;
        qm(j-1,:) = 4*(j-1)/(2*j-3)*chi'.*qm(j,:) - (2*j-1)/(2*j-3)*qm(j+1,:);
    end
    
    % 3. compute error and scale solution
    ratio = q0'./qm(1,:);
    qm = qm.*ratio;

    % 4. compute 1st derivatives
    qmd(1,:) = q0d';
    for i = 1:m
        qmd(i+1,:) = (-(i-1/2)*chi'.*qm(i+1,:) + (i-1/2)*qm(i,:))./(-t'.^2-2*t');
    end
    
    % 5. compute 2nd derivatives
    for i = 1:m+1
        qmdd(i,:) = (-(i-3/2)*(i-1/2)*qm(i,:) + 2*chi'.*qmd(i,:))./(-t'.^2-2*t');
    end
end