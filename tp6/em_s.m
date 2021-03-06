as = load('a.txt');
os = load('o.txt');
us = load('u.txt');

c = {as, os, us};

% leyenda para cada clase
legends = ['a';'o';'u'];

% un color para cada clase
colors = {[1,0,0],[0,1,0],[0,0,1]};

% limites para los graficos
f1_min = 0;
f1_max = 1500;
f2_min = 0;
f2_max = 1800;

% cantidad de clases
K = length(c);

%% mezclar

for i = 1:K
    c{i} = shuffle(c{i});
end

%% me quedo con N muestras y sólo los 2 primeros formantes

f = cell(1, K);
t = cell(1,K);
N = 40;

for k = 1:K
    f{k} = c{k}(1:N, 1:2);
    t{k} = c{k}(N+1:end,1:2);
end

%% valores iniciales
% 
todos = [];

for k = 1:K
    todos = [todos;f{k}];
end

%% esto es para generar muestras random y ver como quedan las secciones

%% seccionar

[xk,u, theta] = seccionar(todos, K);

figure;

for k=1:K
    plot(xk{k}(:,1),xk{k}(:,2),'o', 'color', colors{k});
    hold on;
end

plot(u(1),u(2),'ko','markersize',10, 'markerfacecolor','k');

xlabel('F1 [Hz]');
ylabel('F2 [Hz]');
xlim([f1_min, f1_max])
ylim([f2_min, f2_max])

str = sprintf('%0.1f', theta*180/pi);
title(str);



%% medias y sigmas iniciales para cada clase

% probabilidad de cada clase

p_k = zeros(1,K);
for k=1:K
    p_k(k) = length(xk{k});
end
p_k = p_k/sum(p_k);

u = zeros(K, 2);
for k = 1:K
    u(k, :) = mean(xk{k}, 1);
end

% sigma_k inicial => igual al sigma total de LDA , para todas las clases

sigma = cell(1,K);

for k=1:K
    sigma{k} = calcular_sigma(xk{k},u(k,:));
end

sigma_t = zeros(2,2);

for k=1:K
    sigma_t = sigma_t + p_k(k)*sigma{k};
end

for k=1:K
    sigma{k} = sigma_t;
end

% todas las muestras juntas

xs = [];
for k = 1:K
    xs = [xs;f{k}];
end

NN = length(xs);


%% entrenamiento

p_x = zeros(1,NN); % probabilidad de cada muestra, p(x)
II = 30; % numero de iteraciones
LL = zeros(1,II); % valores de log likelihood para cada iteracion

fig_x = figure;
fig_ll = figure;

for iter = 1:II
    
    gamma_k = zeros(K,NN);
    for i=1:NN
        for k=1:K
            gamma_k(k,i) = mvnpdf(xs(i,:),u(k,:),sigma{k})*p_k(k);
        end
        p_x(i) = sum(gamma_k(:,i));
        gamma_k(:,i) = gamma_k(:,i)/p_x(i);
    end
    
    % pause
    
    % actualizar medias
    
    u = zeros(K,2);
    
    for k=1:K
        u(k,:) = gamma_k(k,:)*xs/sum(gamma_k(k,:));
    end
    
    % actualizar sigmas
    
    sigma = cell(1,K);
    
    for k=1:K
        sigma{k} = zeros(2,2);
        for i=1:NN
            sigma{k} = sigma{k} + gamma_k(k,i)*(xs(i,:)-u(k,:)).'*(xs(i,:)-u(k,:));
        end
        sigma{k} = sigma{k}/sum(gamma_k(k,:));
    end
    
    % actualizar probabilidades de cada clase
    
    for k=1:K
        p_k(k) = sum(gamma_k(k,:))/NN;
    end
    
    % graficar pesando los colores con los gammas
    
    figure(fig_x);
    for i=1:NN
        col = [0 0 0];
        for k=1:K
            col = col + colors{k}*gamma_k(k,i);
        end
        plot(xs(i,1),xs(i,2),'o','color',col,'markerfacecolor',col);
        hold on;
    end
    
    xlabel('F1 [Hz]');
    ylabel('F2 [Hz]');
    xlim([f1_min, f1_max])
    ylim([f2_min, f2_max])
    title('Primeros 2 formantes');
    
    % graficar log likelihood
    
    figure(fig_ll);
    
    LL(iter) = sum(log(p_x));
    
    plot(1:iter,LL(1:iter),'-o');
    xlim([1 II]);
    title('Log Likelihood');
    xlabel('Iteracion');
    
    fprintf('Proxima iteracion ...\n');
    pause(0.1);
    
end

%% test

xs = [];

% pongo todas las muestras juntas

for k = 1:K
    xs = [xs;t{k}];
end

% ws: clase REAL 1..K a la que pertenece cada muestra

ws = [];

for k=1:K
    ws = [ws , ones(1,length(t{k}))*k];
end

% c: clase PREDICHA 1..K para cada muestra

c = clasificar_em(xs,u,sigma,p_k, K);

%% calcular error como #clasificaciones correctas/#total muestras
fprintf('Error: %0.2f %% \n', sum(ws ~= c)/length(xs)*100);

perm = corregir_etiquetas(ws,c, K);

c = perm(c);

% % arreglo los colores
% colors_n = cell(1,K);
% for k=1:K
%     colors_n{k} = colors{perm(k)};
% end
% 
% legends = legends(perm);


fprintf('Error: %0.2f %% \n', sum(ws ~= c)/length(xs)*100);

%% grafico

figure;

% clasificación correcta
graficar_muestras(xs,ws,'x',legends, colors, f1_min, f1_max, f2_min, f2_max, K);

hold on

% resultados de test
graficar_muestras(xs,c,'o',legends, colors, f1_min, f1_max, f2_min, f2_max, K);

title('Test (o) y correcta (x)');

%% Graficar regiones EM

[colors_c, F1, F2] = colores_em(f1_min, f1_max, f2_min, f2_max,5, u, sigma, p_k, K, colors);

%%

figure;
surf(F1,F2,zeros(size(F1)),colors_c,'EdgeColor','none');
xlabel('F1 [Hz]');
ylabel('F2 [Hz]');
view(2);
xlim([f1_min, f1_max])
ylim([f2_min, f2_max])
title('Regiones');


