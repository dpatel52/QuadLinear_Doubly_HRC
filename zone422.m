function [beta,k,M,phi,defl,defl_cant,ectop,es_T,es_C,yield_type] = zone422(beta_z4,L,b,h,alpha,E,epsilon_cr,beta_1,beta_2,beta_3,eta_1,eta_2,eta_3,xi,omega,eta_c,n,kappa,eta_s,rho_c,rho_t)
x=L/2;%location to find deflection along the beam
x2 = L;

k=zeros(size(beta_z4,1),1);
for i=1:size(beta_z4,1)
  beta = beta_z4(i,1);

k(i,1) =  ((-1 + beta_1)*(beta_1 - 2*beta + 1)*eta_1 - (beta_1 - beta_2)*(beta_1 - 2*beta + beta_2)*eta_2 + xi*omega*(beta + omega)*eta_c - (beta_2 - beta)^2*eta_3 - xi*omega*(beta + omega) - n*(eta_s*rho_t + rho_c)*beta^2 + (kappa*n*eta_s*rho_t - kappa*n*rho_t - 2)*beta + sqrt((((eta_s*beta - kappa*(eta_s - 1))*rho_t + beta*rho_c)^2*n^2 + ((2*eta_s*(eta_c*xi*alpha - eta_3*(alpha - 1))*beta^2 + (((4*omega*(eta_c - 1)*alpha + (-2*kappa - 2*omega)*eta_c + 2*omega)*xi - 4*(alpha - 1)*((eta_1 - eta_2)*beta_1 + (eta_2 - eta_3)*beta_2 - eta_1 + 1))*eta_s + 2*eta_c*xi*kappa)*beta + (2*omega*(alpha*omega - kappa - omega)*(eta_c - 1)*xi + 2*(alpha - 1)*((eta_1 - eta_2)*beta_1^2 + (eta_2 - eta_3)*beta_2^2 - eta_1 + 1))*eta_s + 2*xi*kappa*omega*(eta_c - 1))*rho_t - 2*rho_c*((eta_c*(alpha - 1)*xi - alpha*eta_3)*beta^2 + (2*(alpha - 0.5)*omega*(eta_c - 1)*xi - 2*((eta_1 - eta_2)*beta_1 + (eta_2 - eta_3)*beta_2 - eta_1 + 1)*alpha)*beta + alpha*(omega^2*(eta_c - 1)*xi + (eta_1 - eta_2)*beta_1^2 + (eta_2 - eta_3)*beta_2^2 - eta_1 + 1)))*n - xi*(-eta_c*eta_3*beta^2 - 2*eta_c*((eta_1 - eta_2)*beta_1 + (eta_2 - eta_3)*beta_2 - eta_1 + 1)*beta + omega^2*(eta_c - 1)*xi + ((eta_1 - eta_2)*beta_1^2 + (eta_2 - eta_3)*beta_2^2 - eta_1 + 1)*eta_c))*beta^2) + 1)/((-1 + beta_1)*(beta_1 - 2*beta + 1)*eta_1 - (beta_1 - beta_2)*(beta_1 - 2*beta + beta_2)*eta_2 + xi*(beta + omega)^2*eta_c - (beta_2 - beta)^2*eta_3 + (-2*beta*omega - omega^2)*xi - 2*beta + 1);

end

beta=beta_z4;
M = E * b * (((-2 * xi .* eta_c + 2 * eta_3) .* beta.^3 + ((3 * eta_1 - 3 * eta_2) .* beta_1 + (3 * eta_2 - 3 * eta_3) .* beta_2 - 3 * eta_1 + 3 + (-3 * eta_c + 3) .* omega .* xi) .* beta.^2 + (-eta_1 + eta_2) .* beta_1.^3 + (-eta_2 + eta_3) .* beta_2.^3 + eta_1 - 1 + (eta_c - 1) .* omega.^3 .* xi) .* k.^3 + (((-6 * eta_s .* rho_t - 6 * rho_c) .* n - 6 * eta_3) .* beta.^3 + (6 .* rho_t .* kappa .* (eta_s - 1) .* n + (-9 * eta_1 + 9 * eta_2) .* beta_1 + (-9 * eta_2 + 9 * eta_3) .* beta_2 + 9 * eta_1 - 9 + (3 * eta_c - 3) .* omega .* xi) .* beta.^2 + (3 * eta_1 - 3 * eta_2) .* beta_1.^3 + (3 * eta_2 - 3 * eta_3) .* beta_2.^3 - 3 * eta_1 + 3 + (-3 * eta_c + 3) .* omega.^3 .* xi) .* k.^2 + (((12 * eta_s .* rho_t .* alpha - 12 * rho_c .* (alpha - 1)) .* n + 6 * eta_3) .* beta.^3 + (-6 * rho_t .* kappa .* (alpha + 1) .* (eta_s - 1) .* n + (9 * eta_1 - 9 * eta_2) .* beta_1 + (9 * eta_2 - 9 * eta_3) .* beta_2 - 9 * eta_1 + 9) .* beta.^2 + (-3 * eta_1 + 3 * eta_2) .* beta_1.^3 + (-3 * eta_2 + 3 * eta_3) .* beta_2.^3 + 3 * eta_1 - 3 + (3 * eta_c - 3) .* omega.^3 .* xi) .* k + ((-6 * eta_s .* rho_t .* alpha.^2 - 6 * rho_c .* (alpha - 1).^2) .* n - 2 * eta_3) .* beta.^3 + (6 * rho_t .* alpha .* kappa .* (eta_s - 1) .* n + (-3 * eta_1 + 3 * eta_2) .* beta_1 + (-3 * eta_2 + 3 * eta_3) .* beta_2 + 3 * eta_1 - 3) .* beta.^2 + (eta_1 - eta_2) .* beta_1.^3 + (eta_2 - eta_3) .* beta_2.^3 - eta_1 + 1 + (-eta_c + 1) .* omega.^3 .* xi) .* h.^2 .* epsilon_cr ./ (6 .* (k - 1) .* beta.^2);


phi= beta .* epsilon_cr ./ ((-k + 1) .* h);
defl = beta .* x .* epsilon_cr .* (-x + L) ./ (2 .* (k - 1) .* h);
defl_cant = -x2.^2 .* beta .* epsilon_cr ./ (2 .* (-1 + k) .* h);

es_C = -(k - 1 + alpha) .* beta * epsilon_cr ./ (k - 1);
ectop = k .* beta * epsilon_cr ./ (1 - k);
es_T = (-alpha + k) .* beta * epsilon_cr ./ (k - 1);

        if any(es_C >= (kappa*epsilon_cr))
            yield_type = 'top steel yield';
            return;
        else
            yield_type = 'no yield';
        end
end

