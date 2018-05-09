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

for k = 1:K
    c{k} = shuffle(c{k});
end

%% me quedo con N muestras y sólo los 2 primeros formantes

f = cell(1,K);
t = cell(1,K);

N = 40;

for k = 1:K
    f{k} = c{k}(1:N, 1:2);
    t{k} = c{k}(N+1:end,1:2);
end

%% media

todos = [];

for k = 1:K
    todos = [todos;f{k}];
end

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


%%

u = zeros(K, 2);
for k = 1:K
    u(k, :) = mean(xk{k}, 1);
end

%% L iteraciones

c = zeros(1, K * N); % supongo que hay N muestras para cada clase

xs = zeros(K * N, 2);

% pongo todas las muestras juntas en xs

for k = 1:K
    xs((k-1)*N + 1:(k-1)*N + N, :) = f{k};
end

%%
L = 15; % L : cantidad de iteraciones

f_form = figure; % fig de formantes
f_d = figure; % fig de distorsion

d_iter = zeros(1, L); % distorsion total para cada iteracion

for j = 1:L
    d = zeros(1,K); % distorsion para cada clase
    for i = 1:length(xs)
        dist = zeros(1,K);
        for k=1:K
            dist(k) = (xs(i, :) - u(k,:)) * (xs(i, :) - u(k,:))'; % distancia a la media
        end
        [m,l] = min(dist);
        c(i) = l;
        d(l) = d(l) + m;
    end
    d_iter(j) = sum(d);
    
    x_k = cell(1,K);
    
    fprintf('K-Means Iteracion %i\n',j);
    
    for k=1:K
        x_k{k} = xs(c==k,:); % xs clasificados con la clase k en x_k{k}
        fprintf('Total clase %i: %i\n',k,sum(c==k));
    end
    
    fprintf('\n');
    
    % calculo las nuevas medias
    
    for k = 1:K
        u(k, :) = mean(x_k{k}, 1);
    end
    
    % graficos
    
    figure(f_form);
    clf;
    
    for k=1:K
        plot(u(k,1),u(k,2),'o','color',colors{k},'markersize',10,'markerfacecolor',colors{k});
        hold on;
    end
    
    for k=1:K
        plot(x_k{k}(:, 1), x_k{k}(:, 2), 'o','color',colors{k});
        hold on;
    end
    
    
    xlabel('F1 [Hz]');
    ylabel('F2 [Hz]');
    xlim([f1_min, f1_max])
    ylim([f2_min, f2_max])
    
    legend(legends);
    
    str = sprintf('K-Means Iteracion %i', j);
    title(str);
    
    
    figure(f_d);
    clf;
    plot(1:L, d_iter, 'r');
    title(str);
    ylabel('Distorsion');
    
    pause(0.5);
end

% p(k) probabilidad de la clase k

p = zeros(1,K);

for k=1:K
    p(k) = length(x_k{k});
end

p = p/sum(p);

%

f1 = 0:5:2000;
f2 = 0:5:2000;
[F1,F2] = meshgrid(f1,f2);

%%

z = cell(1,K);

sigma_k = cell(1,K);
det_sigma_k = zeros(1,K);
inv_sigma_k = cell(1,K);
surf_k = [];
for k=1:K
    
    u_k = u(k,:);
    sigma_k{k} = calcular_sigma(x_k{k},u_k);
    
    z{k} = zeros(size(F1));
    
    det_sigma_k(k) = det(sigma_k{k}); % calculo esto acá en vez de en g_k
    inv_sigma_k{k} = inv(sigma_k{k}); % así tarda menos
    
    for i=1:numel(F1)
        x = F1(i);
        y = F2(i);
        z{k}(i) = g_k(x,y,det_sigma_k(k),inv_sigma_k{k},u_k, p(k));
    end
    
    surf_k(k) = figure;
    surf(F1,F2,z{k},'EdgeColor','none');
    xlabel('F1 [Hz]');
    ylabel('F2 [Hz]');
    view(2);
    str = sprintf('g_k, k = %i (%s)',k,legends(k));
    xlim([f1_min, f1_max])
    ylim([f2_min, f2_max])
    title(str);
    
end

% esto es para poner el colormap en la misma escala para todos

zz = cell2mat(z);
zz = reshape(zz,1,numel(zz));
min_g_k = min(zz);
max_g_k = max(zz);

for k=1:K
    figure(surf_k(k));
    caxis([min_g_k, max_g_k]);
    colormap jet
end


clas = zeros(size(F1));

for i=1:numel(F1)
    v = zeros(1,K);
    
    for k=1:K
        v(k) = z{k}(i);
    end
    
    [m,clas(i)] = max(v);
end

figure;
h = surf(F1,F2,clas,'EdgeColor','none');
view(2);

xlim([f1_min, f1_max])
ylim([f2_min, f2_max])
map = reshape(cell2mat(colors),K,3);

colormap(map);

% esto es un truco para que aparezca la leyenda en el surf

hold on;

plot_k = [];
for k=1:K
    plot_k(k) = plot(x_k{k}(1,1),x_k{k}(1,2),'o','color',colors{k}, 'markerfacecolor',colors{k});
    hold on;
end

legend(plot_k,legends);
title('Región para cada clase');

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
    
    % p_x(k) = p(k|x)
    
    p_x = zeros(1,K);
    
    for k=1:K
        p_x(k) = g_k(x(1),x(2),det_sigma_k(k), inv_sigma_k{k},u(k,:),p(k));
    end
    
    [m,k_max] = max(p_x);
    
    c(i) = k_max;
end

%% calcular error como #clasificaciones correctas/#total muestras

fprintf('Error: %0.2f %% \n', sum(ws ~= c)/length(xs)*100);

perm = corregir_etiquetas(ws,c, K);

c = perm(c);

fprintf('Error: %0.2f %% \n', sum(ws ~= c)/length(xs)*100);

%% grafico

figure;

% clasificación correcta
graficar_muestras(xs,ws,'x',legends, colors, f1_min, f1_max, f2_min, f2_max, K);

hold on

% resultados de test
graficar_muestras(xs,c,'o',legends, colors, f1_min, f1_max, f2_min, f2_max, K);

title('Test (o) y correcta (x)');

