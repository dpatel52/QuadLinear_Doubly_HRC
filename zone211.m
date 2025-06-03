function [beta,k,M,phi,defl,defl_cant,ectop,es_T,yield_type] = zone211(beta,L,b,h,alpha,E,epsilon_cr,beta_1,beta_2,beta_3,eta_1,eta_2,eta_3,xi,omega,eta_c,n,kappa,eta_s,rho_c,rho_t)
x=L/2;%location to find deflection along the beam
x2=L;
k=real(sqrt((-2 * n .* (-(-1 + beta).^2 .* eta_1 + beta.^2 .* xi - 2 .* beta + 1) .* (rho_c - rho_t) .* alpha + 2 * (n * rho_t + xi / 2) .* (-1 + beta).^2 .* eta_1 + n .* (n .* rho_c^2 + (2 * n .* rho_t + 2 .* xi) .* rho_c + n .* rho_t^2) .* beta.^2 + (4 * n .* rho_t + 2 .* xi) .* beta - 2 * n .* rho_t - xi) .* beta.^2) + (-eta_1 + (-rho_c - rho_t) .* n) .* beta.^2 + (2 .* eta_1 - 2) .* beta - eta_1 + 1) ./ ((xi - eta_1) .* beta.^2 + (2 .* eta_1 - 2) .* beta - eta_1 + 1);

M= -epsilon_cr * E * b * h^2 * ((((-eta_1 / 3 + xi / 3) .* k.^3 + (eta_1 + (rho_c + rho_t) .* n) .* k.^2 + (-eta_1 + ((2 * alpha - 2) * rho_c - 2 * rho_t * alpha) .* n) .* k + eta_1 / 3 + ((alpha^2 - 2 * alpha + 1) * rho_c + rho_t * alpha^2) .* n) .* beta.^3 + (k - 1).^3 .* (eta_1 - 1) .* beta.^2 / 2 - (k - 1).^3 .* (eta_1 - 1) / 6) ./ (beta.^2 .* (k - 1)));


defl = beta .* x .* epsilon_cr .* (-x + L) ./ (2 .* (k - 1) .* h);
defl_cant = -x2.^2 .* beta .* epsilon_cr ./ (2 .* (-1 + k) .* h);


ectop = k .* beta * epsilon_cr ./ (1 - k);
es_T = (-alpha + k) .* beta * epsilon_cr ./ (k - 1);
phi= beta .* epsilon_cr ./ ((-k + 1) .* h);

    if any(ectop >= (omega * epsilon_cr))
        yield_type = 'concrete top yield';
    elseif any(es_T >= (kappa * epsilon_cr))
        yield_type = 'steel yield';
    else
        yield_type = 'no yield';
    end
% 
% % Find the first index where either condition is met
%     index_ectop = find(ectop > (omega * epsilon_cr), 1);
%     index_es_T = find(es_T > (kappa * epsilon_cr), 1);
% 
%     % Determine the stopping index
%     stop_index = min([index_ectop, index_es_T]);
% 
%     if isempty(stop_index)
%         % No threshold is crossed, use all indices
%         stop_index = length(beta);
%     end
% 
%     % Keep only the valid values up to the stopping index
%     beta = beta(1:stop_index);
%     k = k(1:stop_index);
%     M = M(1:stop_index);
%     phi = phi(1:stop_index);
%     defl = defl(1:stop_index);
%     defl_cant = defl_cant(1:stop_index);
%     ectop = ectop(1:stop_index);
%     es_T = es_T(1:stop_index);


end