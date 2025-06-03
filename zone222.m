function [beta,k,M,phi,defl,defl_cant,ectop,es_T] = zone222(beta,L,b,h,alpha,E,epsilon_cr,beta_1,beta_2,beta_3,eta_1,eta_2,eta_3,xi,omega,eta_c,n,kappa,eta_s,rho_c,rho_t)
x=L/2;%location to find deflection along the beam
x2 = L;
k=real((xi .* omega .* (beta + omega) .* eta_c - (-1 + beta).^2 .* eta_1 - xi .* omega.^2 - beta .* omega .* xi - n .* (eta_s .* rho_t + rho_c) .* beta.^2 + (kappa .* n .* eta_s .* rho_t - kappa .* n .* rho_t - 2) .* beta + sqrt((((beta .* eta_s - kappa .* (eta_s - 1)) .* rho_t + rho_c .* beta).^2 .* n.^2 + ((2 .* (xi .* alpha .* eta_c - (alpha - 1) .* eta_1) .* eta_s .* beta.^2 + ((((4 .* alpha .* omega - 2 .* kappa - 2 .* omega) .* eta_c - 4 .* alpha .* omega + 2 .* omega) .* eta_s + 2 .* kappa .* eta_c) .* xi + 4 .* eta_s .* (eta_1 - 1) .* (alpha - 1)) .* beta + 2 .* omega .* ((alpha .* omega - kappa - omega) .* eta_s + kappa) .* (eta_c - 1) .* xi - 2 .* eta_s .* (eta_1 - 1) .* (alpha - 1)) .* rho_t - 2 .* ((eta_c .* (alpha - 1) .* xi - eta_1 .* alpha) .* beta.^2 + (2 .* omega .* (alpha - 1/2) .* (eta_c - 1) .* xi + 2 .* alpha .* (eta_1 - 1)) .* beta + alpha .* (omega.^2 .* (eta_c - 1) .* xi - eta_1 + 1)) .* rho_c) .* n - xi .* (-eta_1 .* eta_c .* beta.^2 + 2 .* eta_c .* (eta_1 - 1) .* beta + omega.^2 .* (eta_c - 1) .* xi - eta_c .* (eta_1 - 1))) .* beta.^2) + 1) ./ (xi .* (beta + omega).^2 .* eta_c - xi .* omega.^2 - 2 .* beta .* omega .* xi - (-1 + beta).^2 .* eta_1 - 2 .* beta + 1));

M= ((((-2 .* xi .* eta_c + 2 .* eta_1) .* beta.^3 + (-3 .* omega .* (eta_c - 1) .* xi - 3 .* eta_1 + 3) .* beta.^2 + omega.^3 .* (eta_c - 1) .* xi + eta_1 - 1) .* k.^3 + (((-6 .* eta_s .* rho_t - 6 .* rho_c) .* n - 6 .* eta_1) .* beta.^3 + (6 .* rho_t .* kappa .* (eta_s - 1) .* n + 3 .* omega .* (eta_c - 1) .* xi + 9 .* eta_1 - 9) .* beta.^2 - 3 .* omega.^3 .* (eta_c - 1) .* xi - 3 .* eta_1 + 3) .* k.^2 + (((12 .* eta_s .* rho_t .* alpha - 12 .* rho_c .* (alpha - 1)) .* n + 6 .* eta_1) .* beta.^3 + (-6 .* rho_t .* kappa .* (alpha + 1) .* (eta_s - 1) .* n - 9 .* eta_1 + 9) .* beta.^2 + 3 .* omega.^3 .* (eta_c - 1) .* xi + 3 .* eta_1 - 3) .* k + ((-6 .* eta_s .* rho_t .* alpha.^2 - 6 .* rho_c .* (alpha - 1).^2) .* n - 2 .* eta_1) .* beta.^3 + (6 .* rho_t .* alpha .* kappa .* (eta_s - 1) .* n + 3 .* eta_1 - 3) .* beta.^2 - omega.^3 .* (eta_c - 1) .* xi - eta_1 + 1) .* b .* E .* epsilon_cr .* h.^2 ./ (6 .* (k - 1) .* beta.^2));


phi= beta .* epsilon_cr ./ ((-k + 1) .* h);
defl = beta .* x .* epsilon_cr .* (-x + L) ./ (2 .* (k - 1) .* h);
defl_cant = -x2.^2 .* beta .* epsilon_cr ./ (2 .* (-1 + k) .* h);

ectop = k .* beta * epsilon_cr ./ (1 - k);
es_T = (-alpha + k) .* beta * epsilon_cr ./ (k - 1);phi= beta .* epsilon_cr ./ ((-k + 1) .* h);



end