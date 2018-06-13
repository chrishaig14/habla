data = load('data.mat');

hmms = {data.hmm1, data.hmm2, data.hmm3, data.hmm4, data.hmm5, data.hmm6};

%% 

h = figure;
l = figure;
m = figure;

%% generar secuencia de observaciones, x

N = 4; % numero del modelo hmm a usar
C = 25; % minimo numero de veces que tiene que aparecer cada estado en la secuencia generada

[x, stateSeq] = generar_x(hmms{N},C);
numStates = length(hmms{N}.means);
[numPts, dim] = size(x);
size(x)

%% inicializar medias, sigmas y matriz de transiciones

hmm = inicializar(x, numStates);

%% entrenamiento

NITER = 100;
n = length(x); % n que voy a usar para calcular la prob a partir de los alphas y betas

logiter = zeros(1,NITER);

logProb_prev = -Inf;

for iter = 1:NITER

    alpha = obtener_alphas(x,hmm);

    beta = obtener_betas(x,hmm);

    gamma = obtener_gammas(alpha, beta);
    
    logProb = logsum(alpha(2:4,n)+beta(2:4,n));
    
    logFwd = logfwd(x,hmm);
    
    fprintf('logProb: %i\n', logProb);
    fprintf('logFwd:  %i\n', logFwd);

    logiter(iter) = logProb;

    xi = obtener_xis(alpha, beta, x, hmm);

    trans = calcular_trans(gamma, xi);

    % para calcular las medias y varianzas, saco el log de los gammas
    
    exp_gamma = exp(gamma);

    u = calcular_mus(exp_gamma, x);
    
    sigma = calcular_sigmas(exp_gamma, x, u);
    
    % si alguna varianza tiene el rcond muy chico, no la cambio
    
    sigma_prev = hmm.vars;
    
    RCONDLIM = 0.005;
    
    for k=2:5-1
        if rcond(sigma{k}) < RCONDLIM
            fprintf('Mantengo sigma anterior: %i -> %i\n', rcond(sigma{k}), rcond(sigma_prev{k}));
            sigma{k} = sigma_prev{k};
        end
    end
    
    % actualizar modelo
    
    hmm.means = u;
    hmm.vars = sigma;
    hmm.trans = trans;

    % actualizar grafico de medias y varianzas
    figure(h)
    clf
    plot(x(:,1),x(:,2),'ob', 'markerfacecolor','b');
    hold on;
    for i=2:4
        plotgaus(hmm.means{i},hmm.vars{i},[1,0,0]);
        hold on;
    end
    d = load('data.mat');
    for i=2:4
        plotgaus(d.hmm4.means{i},d.hmm4.vars{i},[0,1,0]);
        hold on;
    end
    str = sprintf('medias y varianzas reales (g) y estimadas(r), it = %i', iter);
    title(str);

    % actualizar grafico de log-likelihood
    figure(l)
    clf
    plot(1:iter,logiter(1:iter));
    title('log-likelihood');
    xlabel('iteracion');
    
    
    % actualizar grafico de colores con gamas
    figure(m)
    clf
    graficar_gammas(x, gamma);
    str = sprintf('gammas, it = %i', iter);
    title(str);

    % mostrar nueva matriz de transiciones
    
    trans
    
    pause();
    
%   cortar si el log-likelihood cambio "muy poco"
%     if (logProb - logProb_prev) < 1  && iter > 1
%         fprintf('***************************************************\n');
%         fprintf('logProb: %i ; logProb_prev: %i\n;',logProb, logProb_prev);
%         break;
%     end
    
    logProb_prev = logProb;

end
    
