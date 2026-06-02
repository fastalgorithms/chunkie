function [qm, qmd, qmdd] = qleg_half_miller(t, m)

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
    
    % compute chi
    chi = t+1;
    
    % initialize the seeds for the Legendre functions
    [q0,q1,q0d] = chnk.axissymlap2d.qleg_half(t);
    q1d = (-q0 + chi.*q1)/(2*(chi+1)*t);
    
    % check if we can run the forward recurrence
    run_forward = 0;
    if (chi < 1.005)
        run_forward = 1;
    end
    
    qm = zeros(m+1,1);
    qmd = zeros(m+1,1);
    qmdd = zeros(m+1,1);

    % run the forward reccurence
    if run_forward == 1
        qm(1) = q0;
        qm(2) = q1;
        qmd(1) = q0d;
        qmd(2) = q1d;
        for i = 1:(m-1)
            qm(2+i) = 4*(i-1)/(2*i-1)*chi*qm(1+i) - (2*i-3)/(2*i-1)*qm(i);
            qmd(2+i) = (-(m-1/2)*chi.*qm(2+i) + (m+1/2)*qm(1+i))./(1-chi.^2);
        end

        for i = 1:m+1
            qmdd(i) = (-(i-3/2)*(i-1/2)*qm(i) + 2*chi.*qmd(i))./(1-chi.^2);
        end

        return
    end

    % Millers Algorithm:

    % 1. run the forward reccurence until it has blown up
    % recurrence intialization
    fprev = 1;
    fprevprev = 0;
    
    % NOTE:
    % should maybe wait until second derivatives also blow up
    
    maxiter = 100000; % max number of extra terms
    upbound = 1.0e19; % blow up tolerance
    nterms = m+10; % minimum number of extra terms
    for i = m:maxiter
        f = 4*(i-1)/(2*i-1)*chi*fprev - (2*i-3)/(2*i-1)*fprevprev;
        d = (-(m-1/2)*chi.*fprev + (m+1/2)*fprevprev)./(1-chi.^2);

        if (abs(f) >= upbound) 
            if (abs(d) >= upbound)
                nterms = i+1;
                break
            end
        end

        fprevprev = fprev;
        fprev = f;
    end

    % 2. run the backward reccurence
    % recurrence intialization
    fprev = 1;
    fprevprev = 0;
    
    % run the backward reccurence from nterms to m
    for j = 1:(nterms-m+1)
        i = nterms-j+1   +1;

        f = 4*(i-1)/(2*i-3)*chi*fprev - (2*i-1)/(2*i-3)*fprevprev;

        fprevprev = fprev;
        fprev = f;
    end

    qm(m) = fprev;
    qm(m+1) = fprevprev;

    % run the backward reccurence from m to 1
    for j = 1:(m-1)
        i = m-1-j+1   +1;
        qm(i-1) = 4*(i-1)/(2*i-3)*chi*qm(i) - (2*i-1)/(2*i-3)*qm(i+1);
    end

    % 3. compute error and scale solution
    ratio = q0/qm(1);
    qm = qm.*ratio;

    % NOTE:
    % combine these loops to speed things up ?

    % 4. compute 1st derivatives
    qmd(1) = q0d;
    for i = 1:m
        qmd(i+1) = (-(i-1/2)*chi.*qm(i+1) + (i-1/2)*qm(i))./(1-chi.^2);
    end
    
    % 5. compute 2nd derivatives
    for i = 1:m+1
        qmdd(i) = (-(i-3/2)*(i-1/2)*qm(i) + 2*chi.*qmd(i))./(1-chi.^2);
    end
end
%{
function [mask_blowup, mask_stable, nterms] = get_masks(t, m)
    % run the forward reccurence until it has blown up
    % recurrence intialization
    fprev = ones(size(t,1),size(t,2));
    fprevprev = zeros(size(t,1),size(t,2));
    
    % NOTE:
    % should maybe wait until second derivatives also blow up
    
    maxiter = 100000; % max number of extra terms
    ubound = 1.0e20; % upper bound tolerance
    lbound = 1.0e15; % lower bound tolerance
    nterms = m; % minimum number of extra terms
    for i = m:maxiter
        f = 4*(i-1)/(2*i-1)*chi*fprev - (2*i-3)/(2*i-1)*fprevprev;
        d = (-(m-1/2)*chi.*fprev + (m+1/2)*fprevprev)./(1-chi.^2);

        if (max(abs(f),[],"all") >= ubound)
            if (max(abs(d),[],"all") >= ubound)
                nterms = i+1;
                break
            end
        end

        fprevprev = fprev;
        fprev = f;
    end

    if (nterms < 10)
        nterms = 10;
    end

    mask_blowup = (abs(f) >= lbound) & (abs(d) >= lbound);
    mask_stable = ~mask_blowup;
end
%}