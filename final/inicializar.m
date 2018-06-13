function hmm = inicializar(x, numStates)
%% u y sigma iniciales = general para todas las observaciones

u = mean(x,1);
sigma = calcular_sigma(x,u);
u = u';

hmm.means = cell(1,numStates);
hmm.vars = cell(1,numStates);

for k=2:numStates - 1
    hmm.means{k} = u;
    hmm.vars{k} = sigma;
end

hmm.means{1} = [];
hmm.means{numStates} = [];
hmm.vars{1} = [];
hmm.vars{numStates} = [];

hmm.trans = zeros(numStates, numStates);

hmm.trans(1,2) = 1;
for k = 2:numStates - 1
    hmm.trans(k,k)= 0.5;
    hmm.trans(k,k+1)= 0.5;
end

hmm.trans(numStates, numStates) = 1;

end

