function r = g(u_k, sigma, pi_k, x)
    s = inv(sigma);
    r = -1/2*(u_k*s*u_k'-2*(u_k*s)*x') + log(pi_k);
end