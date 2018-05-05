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


%% K : cantidad de clases

K = length(c);

%% mezclar

for i = 1:K
    c{i} = shuffle(c{i});
end

%% me quedo con 40 muestras y sólo los 2 primeros formantes

f = cell(1, K); % train
t = cell(1,K); % test

n_f = 0; % # train
n_t = 0;% # test

N = 40;

for k = 1:K
    f{k} = c{k}(1:N, 1:2);
    n_f = n_f + length(f{k});
    t{k} = c{k}(N+1:end,1:2);
    n_t = n_t + length(t{k});
end

%% graficar primeros dos formantes

figure;
for k=1:K
    plot(f{k}(:,1),f{k}(:,2),'o','color',colors{k});
    hold on;
end
legend(legends);
    xlim([f1_min, f1_max])
    ylim([f2_min, f2_max])
    
title("Datos entrenamiento");

% figure;
% 
% plot(f_a(:,1),f_a(:,2),'ro');
% hold on;
% 
% plot(f_o(:,1),f_o(:,2),'go');
% hold on;
% 
% plot(f_u(:,1),f_u(:,2),'bo');
% hold on;
% 
% axis equal;
% 
% xlabel('F1 [Hz]');
% ylabel('F2 [Hz]');
% 
% legend('a','o','u');

%% calcular u_k inicial para cada clase

u = zeros(K,2);

for k=1:K
    u(k,:) = mean(f{k},1);
end

%% calcular sigma_k para cada clase y después hacer el promedio para sigma total

sigma_k = cell(1,K);

for k=1:K
    sigma_k{k} = calcular_sigma(f{k},u(k,:));
end

% calculo la probabilidad de cada clase como # muestras total clase / # muestras total

p = zeros(1,K);

total = 0;

for k=1:K
    n_k = length(f{k});
    total = total + n_k;
    p(k) = n_k;
end

p = p*1/total;

sigma = zeros(2,2);

for k=1:K
    sigma = sigma + p(k)*sigma_k{k};
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

c = zeros(1,length(xs));

det_sigma = det(sigma);
inv_sigma = inv(sigma);

for i = 1: length(xs)
    x = xs(i,:);
    
    % p_x(k) = p(k|x)
    
    p_x = zeros(1,K);
    
    for k=1:K
        p_x(k) = g_k(x(1),x(2),det_sigma, inv_sigma,u(k,:),p(k));
    end
    
    [m,k_max] = max(p_x);
    
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

%% grafico de g_k para cada clase

f1 = 0:5:2000;
f2 = 0:5:2000;
[F1,F2] = meshgrid(f1,f2);

z = cell(1,K);
surf_k = [];
for k=1:K
    z{k} = zeros(size(F1));
    for i=1:numel(F1)
        x = F1(i);
        y = F2(i);
        z{k}(i) = g_k(x,y,det_sigma, inv_sigma,u(k,:), p(k));
    end
    
    surf_k(k) = figure;
    surf(F1,F2,z{k},'EdgeColor','none');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    view(2)
        xlim([f1_min, f1_max])
    ylim([f2_min, f2_max])
    
    str = sprintf('g_k para k=%i (%s)',k,legends(k));
    title(str)
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



%% grafico de la región correspondiente a cada clase en el plano f1,f2

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
xlabel('F1 [Hz]');
ylabel('F2 [Hz]');
    xlim([f1_min, f1_max])
    ylim([f2_min, f2_max])
    
view(2)

title('Region para cada clase');

hold on

map = reshape(cell2mat(colors),K,3);

colormap(map);

% esto es un truco para que aparezca la leyenda en el surf

plot_k = [];
for k=1:K
    plot_k(k) = plot(x_k{k}(1,1),x_k{k}(1,2),'o','color',colors{k}, 'markerfacecolor',colors{k});
    hold on;
end

legend(plot_k,legends);


%% graficar elipses de varianza
N = 50;
t = [0:2*pi/(N-1):2*pi];
y1 = sin(t);
y2 = cos(t);

y = [y1;y2];

e = cell(1,K);

x = cell(1,K);

figure;

plot_k = [];

for k=1:K
    e{k} = chol(sigma_k{k});
    x{k} = e{k}'*y;
    for i=1:length(x{k})
        x{k}(:,i) = x{k}(:,i) + u(k,:)';
    end
    plot_k(k) = plot(x{k}(1,:),x{k}(2,:),'color',colors{k});
    hold on;
    plot(f{k}(:,1),f{k}(:,2),'o','color',colors{k});
    hold on;
end

title('Varianza para cada clase')

legend(plot_k,legends);
