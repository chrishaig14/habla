function beta = obtener_betas(x,hmm)

means = hmm.means;
vars = hmm.vars;
hmm.trans(hmm.trans<1e-100) = 1e-100;
logTrans = log(hmm.trans);

numStates = length(means);
nMinOne = numStates - 1;
[numPts,dim] = size(x);

log2pi = log(2*pi);
invSig = cell(1,numStates);
logDetVars2 = zeros(1,numStates);

for k=2:nMinOne
  invSig{k} = inv(vars{k});
  logDetVars2(k) = - 0.5 * log(det(vars{k})) - log2pi; % parte de la probabilidad que queda constante, para cada estado. el logaritmo
end

% incializar beta(k,N)

beta = zeros(numStates, numPts);

beta(2:nMinOne,numPts) = logTrans(2:nMinOne,numStates);

% recursion backward

for n = (numPts-1):-1:1
  for k = 2:nMinOne
    X = x(n+1,:)-means{k}';
    logpdf = - 0.5 * (X * invSig{k}) * X' + logDetVars2(k);
    beta(k,n) = logsum(beta(2:nMinOne,n+1)' + logTrans(k,2:nMinOne) + logpdf);
  end
end