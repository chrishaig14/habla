data = load('data.mat');

hmms = {data.hmm1, data.hmm2, data.hmm3, data.hmm4, data.hmm5, data.hmm6};

%% graficar gaussianas para cada estado. son iguales para todos los hmm

mu = hmms{1}.means;
var = hmms{1}.vars;

figure;
for j = 2:length(mu)-1
    plotgaus(mu{j},var{j});
end
title('gaussianas');

%% graficar secuencias de estados generadas para cada hmm

for i = 1:6
    hmm = hmms{i};
    [x,stateSeq] = genhmm(hmm);
    figure;
    plotseq2(x, stateSeq);
    str = sprintf('hmm %i',i);
    title(str);
    
    figure;
    plotseq(x, stateSeq);
    str = sprintf('hmm %i',i);
    title(str);
end

%% 

hmm = hmms{4};

[x,stateSeq] = genhmm(hmm);

alpha = obtener_alphas(x,hmm);
beta = obtener_betas(x,hmm);

[numPts,dim] = size(x);

% p = zeros(1,numPts);
% 
% for n=1:numPts
%     p(n) = logsum(alpha(:,n)+beta(:,n));
% end
% 
% p

gamma = obtener_gammas(alpha, beta);

xi = obtener_xis(alpha, beta, x, hmm);

%% verificar que gamma(k,n) = sum j=1:K xi(j,k)

% for k=1:3
%     for n=2:numPts
% %         xi(:,k,n)
%         fprintf('logsum:\t %i\n', logsum(xi(:,k,n)));
%         fprintf('gamma:\t %i\n\n', gamma(k,n));
%     end
% end

%% juntar hmm4 con hmm6

hmm4 = hmms{4};
hmm6 = hmms{6};

hmm4_t = hmm4.trans;
hmm6_t = hmm6.trans;

hmm46_t = zeros(8,8);

% en hmm46:
% estado 1: inicial
% estado 2: primero de hmm4
% estado 4: ultimo de hmm4
% estado 5: primero de hmm6
% estado 7: ultimo de hmm6
% estado 8: final

hmm46_t(1,2) = 0.5; % probabilidad de empezar en hmm4
hmm46_t(1,5) = 0.5; % probabilidad de empezar en hmm6
hmm46_t(2:4,2:4) = hmm4_t(2:4,2:4); % dentro de hmm4 o hmm6 las probabilidades son las mismas
hmm46_t(5:7,5:7) = hmm6_t(2:4,2:4);

hmm46_t(4,[2,5,8]) = hmm4_t(4,5)/3; % reparto la probabilidad entre ir a hmm4, hmm5 o al final
hmm46_t(7,[2,5,8]) = hmm6_t(4,5)/3; 

hmm46_t(8,8) = 1;

hmm46_m = hmm4.means;
hmm46_m{5} = hmm6.means{2};
hmm46_m{6} = hmm6.means{3};
hmm46_m{7} = hmm6.means{4};
hmm46_m{8} = []; % para el estado final

hmm46_v = hmm4.vars;
hmm46_v{5} = hmm6.vars{2};
hmm46_v{6} = hmm6.vars{3};
hmm46_v{7} = hmm6.vars{4};
hmm46_v{8} = []; % para el estado final

hmm46 = struct();
hmm46.trans = hmm46_t;
hmm46.vars = hmm46_v;
hmm46.means = hmm46_m;

%% generar secuencia en hmm46
    
    C = 10; % cada estado tiene que aparecer mínimo C veces

    [x, stateSeq_gen] = generar_x(hmm46, C);
    
    seq_gen = obtener_sec_modelos(stateSeq_gen);
    
    [stateSeq_vit, logProb_vit] = logvit(x,hmm46);
    
    seq_vit = obtener_sec_modelos(stateSeq_vit);
    
    fprintf("Secuencia de estados generada: \n");
    stateSeq_gen
    fprintf("Secuencia de estados Viterbi: \n");
    stateSeq_vit
    
    fprintf('Las secuencias de estados coinciden en %i/%i\n',sum(stateSeq_gen==stateSeq_vit), length(stateSeq_vit));
    
    fprintf("Secuencia de modelos generada: \n");
    seq_gen
    fprintf("Secuencia de modelos Viterbi: \n");
    seq_vit
    
    alpha = obtener_alphas(x,hmm46);
    beta = obtener_betas(x,hmm46);
    
    n = length(x);
    
    logProb_gen = logsum(alpha(2:numStates-1,n)+beta(2:numStates-1,n));
    
    fprintf("logProb real = %i\n", logProb_gen);
    fprintf("logProb Viterbi = %i\n", logProb_vit);
    