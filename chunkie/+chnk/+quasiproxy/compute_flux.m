function [flux_up, flux_down] = compute_flux(au, ad, k1x, k1, kN, d, BR)
% - Flux error estimator, from Cho and Barnett's solver


% flux_up = 0;
% flux_down = 0;

flux_up = zeros(2*BR+1,1);
flux_down = zeros(2*BR+1,1);
for J=-BR:BR
    kappa_n  = k1x+2*pi*J/(d);
    ku_n     = sqrt(k1^2-kappa_n^2);
    kd_n     = sqrt(kN^2-kappa_n^2);
    
    flux_up(J+BR+1) = ku_n*(abs(au(J+BR+1))^2);
    flux_down(J+BR+1)=kd_n*(abs(ad(J+BR+1))^2);
%     flux_up  =  flux_up+ku_n*(abs(au(J+BR+1))^2);
%     flux_down = flux_down+kd_n*(abs(ad(J+BR+1))^2);
end

