k = 5;
[alpha_exact,beta_exact] = max_order_explicit_coefficients(k);
rho_exact = fliplr(alpha_exact);
sigma_exact = fliplr(beta_exact);

% controlla se il metodo è 0-stabile
zero_stable = is_von_neumann(rho_exact);
