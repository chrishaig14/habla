function alpha = obtener_alphas(x,hmm)

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

% Initialize the alpha vector for the emitting states

alpha = zeros(numStates, numPts);

for k=2:nMinOne
  X = x(1,:)-means{k}';
  alpha(k,1) = logTrans(1,k) ...
      - 0.5 * (X * invSig{k}) * X' + logDetVars2(k);
end

% Do the forward recursion
for n = 2:numPts
  for k = 2:nMinOne
    X = x(n,:)-means{k}';
    logpdf = - 0.5 * (X * invSig{k}) * X' + logDetVars2(k);
    alpha(k,n) = logsum(alpha(2:nMinOne,n-1) + logTrans(2:nMinOne,k) ) + logpdf;
  end
end