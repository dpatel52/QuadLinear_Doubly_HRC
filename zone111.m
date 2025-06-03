function [beta,k,M,phi,defl,defl_cant] = zone111(beta,L,b,h,alpha,E,epsilon_cr,beta_1,beta_2,beta_3,eta_1,eta_2,eta_3,xi,omega,eta_c,n,kappa,eta_s,rho_c,rho_t)
x=L/2;%location to find deflection along the beam
x2=L;
k=ones(size(beta))*(sqrt((rho_t + rho_c)^2*n^2 + (((-2*xi + 2)*alpha + 2*xi)*rho_c + 2*rho_t*(1 + (xi - 1)*alpha))*n + xi) - 1 + (-rho_c - rho_t)*n)/(xi - 1);
M= -3 * E * b .* beta .* epsilon_cr .* ((xi/3 - 1/3) .* k.^3 + (1 + (rho_t + rho_c) .* n) .* k.^2 + (-1 + ((2*alpha - 2) .* rho_c - 2*rho_t .* alpha) .* n) .* k + 1/3 + ((alpha - 1).^2 .* rho_c + rho_t .* alpha.^2) .* n) .* h.^2 ./ (3 .* k - 3);

phi=beta .* epsilon_cr ./ ((-k + 1) .* h);
defl = beta .* x .* epsilon_cr .* (-x + L) ./ (2 .* (k - 1) .* h);
defl_cant = -x2.^2 .* beta .* epsilon_cr ./ (2 .* (-1 + k) .* h);
end