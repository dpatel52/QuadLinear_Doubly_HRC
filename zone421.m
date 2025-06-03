function [beta,k,M,phi,defl,defl_cant,ectop,es_T,yield_type] = zone421(beta,L,b,h,alpha,E,epsilon_cr,beta_1,beta_2,beta_3,eta_1,eta_2,eta_3,xi,omega,eta_c,n,kappa,eta_s,rho_c,rho_t)
x=L/2;%location to find deflection along the beam
x2 = L;
k = real((sqrt((-2 * (-1 + beta_1) .* (beta_1 - 2 .* beta + 1) .* (xi .* eta_c / 2 + ((rho_c - rho_t) .* alpha + rho_t) .* n) .* eta_1 + 2 * (xi .* eta_c / 2 + ((rho_c - rho_t) .* alpha + rho_t) .* n) .* (beta_1 - 2 .* beta + beta_2) .* (beta_1 - beta_2) .* eta_2 - 2 * xi .* (-(beta_2 - beta).^2 .* eta_3 / 2 + xi .* omega.^2 / 2 + n .* (beta + omega).^2 .* (rho_c - rho_t) .* alpha + n .* rho_t .* omega.^2 - n .* beta .* (rho_c - rho_t) .* omega - rho_c .* n .* beta.^2 - beta + 1 / 2) .* eta_c + 2 * ((rho_c - rho_t) .* alpha + rho_t) .* n .* (beta_2 - beta).^2 .* eta_3 + omega.^2 .* xi.^2 + 4 * n .* omega .* ((beta + omega / 2) .* (rho_c - rho_t) .* alpha + rho_t .* omega / 2 - beta .* (rho_c - rho_t) / 2) .* xi + n .* (4 * (beta - 1 / 2) .* (rho_c - rho_t) .* alpha + n .* (rho_c + rho_t).^2 .* beta.^2 + 4 .* beta .* rho_t - 2 .* rho_t)) .* beta.^2) + (-eta_3 + (-rho_c - rho_t) .* n) .* beta.^2 + (omega .* (eta_c - 1) .* xi + (-2 .* eta_1 + 2 .* eta_2) .* beta_1 + (-2 .* eta_2 + 2 .* eta_3) .* beta_2 + 2 .* eta_1 - 2) .* beta + omega.^2 .* (eta_c - 1) .* xi + (eta_1 - eta_2) .* beta_1.^2 + (eta_2 - eta_3) .* beta_2.^2 - eta_1 + 1) ./ ((xi .* eta_c - eta_3) .* beta.^2 + (2 .* omega .* (eta_c - 1) .* xi + (-2 .* eta_1 + 2 .* eta_2) .* beta_1 + (-2 .* eta_2 + 2 .* eta_3) .* beta_2 + 2 .* eta_1 - 2) .* beta + omega.^2 .* (eta_c - 1) .* xi + (eta_1 - eta_2) .* beta_1.^2 + (eta_2 - eta_3) .* beta_2.^2 - eta_1 + 1));

M = -E * epsilon_cr * (((2 * xi * eta_c - 2 * eta_3) .* beta.^3 + ((-3 * eta_1 + 3 * eta_2) .* beta_1 + (-3 * eta_2 + 3 * eta_3) .* beta_2 + 3 * eta_1 - 3 + (3 * eta_c - 3) .* omega .* xi) .* beta.^2 + (eta_1 - eta_2) .* beta_1.^3 + (eta_2 - eta_3) .* beta_2.^3 - eta_1 + 1 + (-eta_c + 1) .* omega.^3 .* xi) .* k.^3 + ((6 * eta_3 + (6 * rho_c + 6 * rho_t) .* n) .* beta.^3 + ((9 * eta_1 - 9 * eta_2) .* beta_1 + (9 * eta_2 - 9 * eta_3) .* beta_2 - 9 * eta_1 + 9 + (-3 * eta_c + 3) .* omega .* xi) .* beta.^2 + (-3 * eta_1 + 3 * eta_2) .* beta_1.^3 + (-3 * eta_2 + 3 * eta_3) .* beta_2.^3 + 3 * eta_1 - 3 + (3 * eta_c - 3) .* omega.^3 .* xi) .* k.^2 + ((-6 * eta_3 + ((12 * alpha - 12) .* rho_c - 12 .* rho_t .* alpha) .* n) .* beta.^3 + ((-9 * eta_1 + 9 * eta_2) .* beta_1 + (-9 * eta_2 + 9 * eta_3) .* beta_2 + 9 * eta_1 - 9) .* beta.^2 + (3 * eta_1 - 3 * eta_2) .* beta_1.^3 + (3 * eta_2 - 3 * eta_3) .* beta_2.^3 - 3 * eta_1 + 3 + (-3 * eta_c + 3) .* omega.^3 .* xi) .* k + (2 * eta_3 + ((6 * alpha.^2 - 12 * alpha + 6) .* rho_c + 6 * rho_t .* alpha.^2) .* n) .* beta.^3 + ((3 * eta_1 - 3 * eta_2) .* beta_1 + (3 * eta_2 - 3 * eta_3) .* beta_2 - 3 * eta_1 + 3) .* beta.^2 + (-eta_1 + eta_2) .* beta_1.^3 + (-eta_2 + eta_3) .* beta_2.^3 + eta_1 - 1 + (eta_c - 1) .* omega.^3 .* xi) .* b .* h.^2 ./ (6 * (k - 1) .* beta.^2);


phi= beta .* epsilon_cr ./ ((-k + 1) .* h);
defl = beta .* x .* epsilon_cr .* (-x + L) ./ (2 .* (k - 1) .* h);
defl_cant = -x2.^2 .* beta .* epsilon_cr ./ (2 .* (-1 + k) .* h);

ectop = k .* beta * epsilon_cr ./ (1 - k);
es_T = (-alpha + k) .* beta * epsilon_cr ./ (k - 1);

        if any(es_T >= (kappa*epsilon_cr))
            yield_type = 'steel yield';
        else
            yield_type = 'no yield';
        end

% 
% % Logical indexing to filter results based on the threshold values
% valid_indices = (es_T <= (kappa*epsilon_cr));
% 
% beta=beta(valid_indices);
% k = k(valid_indices);
% M = M(valid_indices);
% phi = phi(valid_indices);
% defl = defl(valid_indices);
% defl_cant = defl_cant(valid_indices);
% ectop = ectop(valid_indices);
% es_T = es_T(valid_indices);



end