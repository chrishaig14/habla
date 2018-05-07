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

%p_k inicial => equiprobable

p_k = ones(1,K)*1/K;

%% media

% calculo la media general

todos = [];

for k = 1:K
    todos = [todos;f{k}];
end

u = mean(todos,1);

figure;
plot(todos(:,1),todos(:,2),'o');
hold on;
plot(u(1),u(2),'*');


%%

%figure;


% K = 3;
% colors = {'r','g','b','k','m','y','c'};
%
% todos = zeros(1000,2);
%
% for i=1:length(todos)
%     todos(i,1) = rand()*(f1_max-f1_min)+f1_min;
%     todos(i,2) = rand()*(f2_max-f2_min)+f2_min;
% end

sec = figure;
%%

theta = rand()*2*pi/K; % medido desde la horizontal en u

theta


figure(sec);
plot(u(1),u(2),'*');
hold on;


thetas = [theta:2*pi/K:theta + 2*pi/K*(K-1)];

thetas = wrapTo2Pi(thetas);

xk = cell(1,K);

for i=1:length(todos)
    x_o = todos(i,:); % un punto
    
    % lo llevo hasta u
    
    x = x_o - u;
    
    angulo2pi = wrapTo2Pi(atan2(x(2),x(1))); % angulo entre x y u, desde la horizontal
    
    for k = 1:K-1
        angulo = angulo2pi - thetas(k);
        fin = thetas(k+1) - thetas(k);
        angulo = wrapTo2Pi(angulo);
        fin = wrapTo2Pi(fin);
        if angulo < fin
            xk{k} = [xk{k};x_o];
            %plot(x_o(1),x_o(2),'o', 'color', colors{k});
        end
    end
    angulo = angulo2pi - thetas(K);
    fin = thetas(1) - thetas(K);
    angulo = wrapTo2Pi(angulo);
    fin = wrapTo2Pi(fin);
    if angulo < fin
        xk{K} = [xk{K};x_o];
        %         plot(x_o(1),x_o(2),'o', 'color', colors{K});
    end
    
    hold on;
end
for k=1:K
    plot(xk{k}(:,1),xk{k}(:,2),'o', 'color', colors{k});
    hold on;
end

xlabel('F1 [Hz]');
ylabel('F2 [Hz]');
str = sprintf("%0.1f", 360/K);
title(str);
xlim([f1_min, f1_max])
ylim([f2_min, f2_max])



%%

%M = 1; % con M más grande es muy "fácil"
p_k = zeros(1,K);
for k=1:K
    p_k(k) = length(xk{k});
end
p_k = p_k/sum(p_k);

u = zeros(K, 2);
for k = 1:K
    %u(k, :) = mean(f{k}(1:M, :), 1);
    u(k, :) = mean(xk{k}, 1);
    %u(k,:)
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
