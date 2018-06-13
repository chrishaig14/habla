function [] = verificar(alpha, beta, gamma, xi)

    [K, numPts] = size(alpha);

% verificar alphas y betas

    p = zeros(1,numPts);

    for n=1:numPts
        p(n) = logsum(alpha(:,n)+beta(:,n));
    end

    p

% verificar gammas y xis

    r=zeros(3,numPts-1);

    for n=1:numPts-1
        for j=1:3
            r(j,n) = gamma(j,n)- logsum(xi(j,:,n));
%             r(j,n) = exp(gamma(j,n))/sum(exp(xi(j,:,n)));
        end
    end

    r

% verificar gammas y xis

    r=zeros(3,numPts-1);

    for n=1:numPts-1
        for k=1:3
                        r(k,n) = gamma(k,n+1)-logsum(xi(:,k,n));
%             r(k,n) = exp(gamma(k,n+1))/sum(exp(xi(:,k,n)));
        end
    end

    r
%
end