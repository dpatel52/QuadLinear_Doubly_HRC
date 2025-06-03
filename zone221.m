function [beta,k,M,phi,defl,defl_cant,ectop,es_T,yield_type] = zone221(beta,L,b,h,alpha,E,epsilon_cr,beta_1,beta_2,beta_3,eta_1,eta_2,eta_3,xi,omega,eta_c,n,kappa,eta_s,rho_c,rho_t)
x=L/2;%location to find deflection along the beam
x2=L;

k= real((-(-1 + beta).^2 .* eta_1 + xi .* omega .* (beta + omega) .* eta_c - n .* (rho_c + rho_t) .* beta.^2 + (-omega .* xi - 2) .* beta - xi .* omega.^2 + sqrt(beta.^2 .* (beta.^2 .* (rho_c + rho_t).^2 .* n.^2 + ((-2 .* eta_c .* (-rho_t .* alpha + rho_c .* (alpha - 1)) .* xi + 2 .* eta_1 .* ((-alpha + 1) .* rho_t + rho_c .* alpha)) .* beta.^2 + (-4 .* (alpha - 1/2) .* omega .* (eta_c - 1) .* (rho_c - rho_t) .* xi - 4 .* (eta_1 - 1) .* ((-alpha + 1) .* rho_t + rho_c .* alpha)) .* beta - 2 .* ((-alpha + 1) .* rho_t + rho_c .* alpha) .* (omega.^2 .* (eta_c - 1) .* xi - eta_1 + 1)) .* n - xi .* (-eta_1 .* eta_c .* beta.^2 + 2 .* eta_c .* (eta_1 - 1) .* beta + omega.^2 .* (eta_c - 1) .* xi - eta_c .* (eta_1 - 1)))) + 1) ./ (xi .* (beta + omega).^2 .* eta_c - (-1 + beta).^2 .* eta_1 - 2 .* beta .* omega .* xi - xi .* omega.^2 - 2 .* beta + 1));

M=(b * ((((-xi .* eta_c + eta_1) .* beta.^3 + (-(3 .* omega .* (eta_c - 1) .* xi) ./ 2 + 3 ./ 2 - (3 .* eta_1) ./ 2) .* beta.^2 + omega.^3 .* (eta_c - 1) .* xi ./ 2 - 1 ./ 2 + eta_1 ./ 2) .* k.^3 + ((-3 .* eta_1 + (-3 .* rho_c - 3 .* rho_t) .* n) .* beta.^3 + ((3 .* omega .* (eta_c - 1) .* xi) ./ 2 - 9 ./ 2 + (9 .* eta_1) ./ 2) .* beta.^2 - (3 .* omega.^3 .* (eta_c - 1) .* xi) ./ 2 + 3 ./ 2 - (3 .* eta_1) ./ 2) .* k.^2 + ((3 .* eta_1 + ((-6 .* alpha + 6) .* rho_c + 6 .* rho_t .* alpha) .* n) .* beta.^3 + (9 ./ 2 - (9 .* eta_1) ./ 2) .* beta.^2 + (3 .* omega.^3 .* (eta_c - 1) .* xi) ./ 2 - 3 ./ 2 + (3 .* eta_1) ./ 2) .* k + (-eta_1 + ((-3 .* alpha.^2 + 6 .* alpha - 3) .* rho_c - 3 .* rho_t .* alpha.^2) .* n) .* beta.^3 + (-3 ./ 2 + (3 .* eta_1) ./ 2) .* beta.^2 - omega.^3 .* (eta_c - 1) .* xi ./ 2 + 1 ./ 2 - eta_1 ./ 2) .* E .* epsilon_cr .* h.^2 ./ (3 .* (k - 1) .* beta.^2)));


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

