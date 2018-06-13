function xi = obtener_xis(alpha, beta, x, hmm)

means = hmm.means;
vars = hmm.vars;
hmm.trans(hmm.trans<1e-100) = 1e-100;
logTrans = log(hmm.trans);

% logTrans

log2pi = log(2*pi);

numStates = length(means);
nMinOne = numStates - 1;

invSig = cell(1,numStates);
logDetVars2 = zeros(1,numStates);

for k=2:nMinOne
  invSig{k} = inv(vars{k});
  logDetVars2(k) = - 0.5 * log(det(vars{k})) - log2pi;
end

numPts = length(x);

xi = zeros(numStates,numStates, numPts); % LOG de xi!

for n = 2:numPts
    for j = 2:nMinOne
        for k = 2:nMinOne
            X = x(n,:)-means{k}';
            logpdf = - 0.5 * (X * invSig{k}) * X' + logDetVars2(k);
%             fprintf('logpdf x(%i) = %i, k = %i\n',n,logpdf,k);
            xi(j,k,n) = alpha(j,n-1) +  beta(k,n) + logTrans(j,k) + logpdf;
        end
    end
    xi(2:nMinOne,2:nMinOne,n) = xi(2:nMinOne,2:nMinOne,n)-logsum(alpha(2:nMinOne,n)+beta(2:nMinOne,n));
end