as = load('a.txt');
os = load('o.txt');
us = load('u.txt');

c = {as, os, us};

legends = ['a';'o';'u'];
colors = {[1,0,0],[0,1,0],[0,0,1]};

f1_min = 0;
f1_max = 1500;
f2_min = 0;
f2_max = 1800;


K = length(c); % cantidad de clases

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

p_k = zeros(1,K);
for k=1:K
    p_k(k) = length(f{k});
end
p_k = p_k/sum(p_k);

% u_k inicial => tomo M de cada clase y calculo la media

u = zeros(K,2);

M = 5;

for i = 1:K
    u(i, :) = mean(f{i}(1:M, :), 1);
end

% sigma_k inicial => igual al sigma total de LDA , para todas las clases

sigma = cell(1,K);

for k=1:K
    sigma{k} = calcular_sigma(f{k},u(k,:));
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


%%

p_x = zeros(1,NN);
II = 50;
LL = zeros(1,II);

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
    
    %pause
    u = zeros(K,2);
    
    for k=1:K
        u(k,:) = gamma_k(k,:)*xs/sum(gamma_k(k,:));
    end
    
    sigma = cell(1,K);
    
    for k=1:K
        sigma{k} = zeros(2,2);
        for i=1:NN
            sigma{k} = sigma{k} + gamma_k(k,i)*(xs(i,:)-u(k,:)).'*(xs(i,:)-u(k,:));
        end
        sigma{k} = sigma{k}/sum(gamma_k(k,:));
    end
    
    for k=1:K
        p_k(k) = sum(gamma_k(k,:))/NN;
    end
    
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
    
    figure(fig_ll);
    
    LL(iter) = sum(log(p_x));
    
    plot(1:iter,LL(1:iter),'-o');
    xlim([1 II]);
    title('Log Likelihood');
    xlabel('Iteracion');
    
    fprintf('Proxima iteracion ...\n');
    pause();
    
end

%% test

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

%% graficos

% x_k{k} muestras clasificadas para la clase k

x_k = cell(1,K);

for k=1:K
    x_k{k} = xs(c==k,:);
end

figure;

for k=1:K
    plot(x_k{k}(:,1),x_k{k}(:,2),'o','color',colors{k});
    hold on;
end
xlim([f1_min, f1_max])
ylim([f2_min, f2_max])

legend(legends);

title('Resultados test');


%% calcular error como #clasificaciones correctas/#total muestras
fprintf('Error: %0.2f %% \n', sum(ws ~= c)/length(xs)*100);
