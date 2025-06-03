function [beta,k,M,phi,defl,defl_cant,ectop,es_T,yield_type] = zone311(beta,L,b,h,alpha,E,epsilon_cr,beta_1,beta_2,beta_3,eta_1,eta_2,eta_3,xi,omega,eta_c,n,kappa,eta_s,rho_c,rho_t)
x=L/2;%location to find deflection along the beam
x2 = L;

k=real((sqrt(beta.^2 .* (4 .* (n .* (rho_c - rho_t) .* alpha + n .* rho_t + xi / 2) .* (-1 + beta_1) .* (-beta_1 / 2 + beta - 1 / 2) .* eta_1 + 2 .* (beta - beta_1).^2 .* (n .* (rho_c - rho_t) .* alpha + n .* rho_t + xi / 2) .* eta_2 - 2 .* n .* (rho_c - rho_t) .* (beta.^2 .* xi - 2 .* beta + 1) .* alpha + (n .* rho_c.^2 + (2 .* n .* rho_t + 2 .* xi) .* rho_c + n .* rho_t.^2) .* n .* beta.^2 + (4 .* n .* rho_t + 2 .* xi) .* beta - 2 .* n .* rho_t - xi)) + (-eta_2 + (-rho_c - rho_t) .* n) .* beta.^2 + ((-2 .* eta_1 + 2 .* eta_2) .* beta_1 + 2 .* eta_1 - 2) .* beta + (-eta_2 + eta_1) .* beta_1.^2 - eta_1 + 1) ./ ((xi - eta_2) .* beta.^2 + ((-2 .* eta_1 + 2 .* eta_2) .* beta_1 + 2 .* eta_1 - 2) .* beta + (-eta_2 + eta_1) .* beta(1).^2 - eta_1 + 1));

M= h^2 * b * (((2 * eta_2 - 2 * xi) .* beta.^3 + ((-3 * eta_2 + 3 * eta_1) .* beta_1 - 3 * eta_1 + 3) .* beta.^2 + (eta_2 - eta_1) .* beta_1.^3 + eta_1 - 1) .* k.^3 + ((-6 * eta_2 + (-6 * rho_c - 6 * rho_t) .* n) .* beta.^3 + ((9 * eta_2 - 9 * eta_1) .* beta_1 + 9 * eta_1 - 9) .* beta.^2 + (-3 * eta_2 + 3 * eta_1) .* beta_1.^3 - 3 * eta_1 + 3) .* k.^2 + ((6 * eta_2 + ((-12 * alpha + 12) .* rho_c + 12 .* rho_t .* alpha) .* n) .* beta.^3 + ((-9 * eta_2 + 9 * eta_1) .* beta_1 - 9 * eta_1 + 9) .* beta.^2 + (3 * eta_2 - 3 * eta_1) .* beta_1.^3 + 3 * eta_1 - 3) .* k + (-2 * eta_2 + ((-6 * alpha.^2 + 12 * alpha - 6) .* rho_c - 6 .* rho_t .* alpha.^2) .* n) .* beta.^3 + ((3 * eta_2 - 3 * eta_1) .* beta_1 + 3 * eta_1 - 3) .* beta.^2 + (-eta_2 + eta_1) .* beta_1.^3 - eta_1 + 1) .* epsilon_cr .* E ./ (6 * beta.^2 .* (k - 1));


ectop = k .* beta * epsilon_cr ./ (1 - k);
es_T = (-alpha + k) .* beta * epsilon_cr ./ (k - 1);
phi= beta .* epsilon_cr ./ ((-k + 1) .* h);
defl = beta .* x .* epsilon_cr .* (-x + L) ./ (2 .* (k - 1) .* h);
defl_cant = -x2.^2 .* beta .* epsilon_cr ./ (2 .* (-1 + k) .* h);


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