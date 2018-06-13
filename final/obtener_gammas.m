function gamma = obtener_gammas(alpha, beta)

[numStates, numPts] = size(alpha);
nMinOne = numStates - 1;
gamma = zeros(numStates, numPts); % LOG de gamma!

for n = 1:numPts
    for k = 2:nMinOne
        gamma(k,n) = alpha(k,n)+beta(k,n);
    end
    gamma(2:nMinOne,n) = gamma(2:nMinOne,n)-logsum(alpha(2:nMinOne,n)+beta(2:nMinOne,n));
end