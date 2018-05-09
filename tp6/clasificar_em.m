function [ c ] = clasificar_em( xs, u, sigma, p_k, K )
    c = zeros(1,length(xs));

    for i = 1: length(xs)
        x = xs(i,:);
        gamma_k = zeros(K);
        for k=1:K
            gamma_k(k) = mvnpdf(xs(i,:),u(k,:),sigma{k})*p_k(k);
        end
        p_x = sum(gamma_k);
        gamma_k = gamma_k/p_x;

        [m,k_max] = max(gamma_k);

        c(i) = k_max;
    end
end

