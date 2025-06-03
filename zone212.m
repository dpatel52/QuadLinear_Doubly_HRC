function [beta,k,M,phi,defl,defl_cant,ectop,es_T,yield_type] = zone212(beta,L,b,h,alpha,E,epsilon_cr,beta_1,beta_2,beta_3,eta_1,eta_2,eta_3,xi,omega,eta_c,n,kappa,eta_s,rho_c,rho_t)
x=L/2;%location to find deflection along the beam
x2=L;

k=real(-(-1 + beta).^2 .* eta_1 - n .* (eta_s .* rho_t + rho_c) .* beta.^2 + (kappa .* n .* eta_s .* rho_t - kappa .* n .* rho_t - 2) .* beta + sqrt(beta.^2 .* (((eta_s .* rho_t + rho_c) .* beta - rho_t .* kappa .* (eta_s - 1)).^2 .* n.^2 + ((2 .* ((-alpha + 1) .* eta_1 + xi .* alpha) .* eta_s .* rho_t - 2 .* rho_c .* (-alpha .* eta_1 + xi .* (alpha - 1))) .* beta.^2 + ((((4 .* alpha - 4) .* eta_1 - 2 .* xi .* kappa - 4 .* alpha + 4) .* eta_s + 2 .* xi .* kappa) .* rho_t - 4 .* alpha .* rho_c .* (eta_1 - 1)) .* beta - 2 .* ((alpha - 1) .* eta_s .* rho_t - rho_c .* alpha) .* (eta_1 - 1)) .* n + xi .* (beta.^2 .* eta_1 + (-2 .* eta_1 + 2) .* beta + eta_1 - 1))) + 1) ./ (-(-1 + beta).^2 .* eta_1 + beta.^2 .* xi - 2 .* beta + 1);

M= h^2 * E * epsilon_cr * (((eta_1 - xi) .* k.^3 + ((-3 * eta_s * rho_t - 3 * rho_c) .* n - 3 * eta_1) .* k.^2 + ((6 * eta_s * rho_t * alpha - 6 * rho_c * (alpha - 1)) .* n + 3 * eta_1) .* k + (-3 * eta_s * rho_t * alpha.^2 - 3 * rho_c * (alpha - 1).^2) .* n - eta_1) .* beta.^3 - 3 * ((eta_1 - 1) .* k.^2 + (-2 * rho_t * kappa * (eta_s - 1) .* n - 2 * eta_1 + 2) .* k + 2 * rho_t * alpha * kappa * (eta_s - 1) .* n + eta_1 - 1) .* (k - 1) .* beta.^2 / 2 + (k - 1).^3 .* (eta_1 - 1) / 2) .* b ./ (3 * beta.^2 .* (k - 1));


phi= beta .* epsilon_cr ./ ((-k + 1) .* h);
defl = beta .* x .* epsilon_cr .* (-x + L) ./ (2 .* (k - 1) .* h);
defl_cant = -x2.^2 .* beta .* epsilon_cr ./ (2 .* (-1 + k) .* h);

ectop = k .* beta * epsilon_cr ./ (1 - k);
es_T = (-alpha + k) .* beta * epsilon_cr ./ (k - 1);

        if any(ectop >= (omega * epsilon_cr))
            yield_type = 'concrete top yield';
        else
            yield_type = 'no yield';
        end

    % % Logical indexing to filter results based on the threshold values
    % valid_indices = (ectop <= (omega*epsilon_cr));
    % 
    % beta=beta(valid_indices);
    % k = k(valid_indices);
    % M = M(valid_indices);
    % phi = phi(valid_indices);
    % defl = defl(valid_indices);
    % defl_cant = defl_cant(valid_indices);
    % ectop = ectop(valid_indices);
    % es_T = es_T(valid_indices);
    % 
end