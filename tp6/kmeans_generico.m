as = load('a.txt');
os = load('o.txt');
us = load('u.txt');

c = {as, os, us};

col = ['r', 'g', 'b']; % colores para los graficos, uno para cada clase

K = length(c); % cantidad de clases

%% mezclar

for i = 1:K
    c{i} = shuffle(c{i});
end

%% me quedo con N muestras y sólo los 2 primeros formantes

f = cell(1, K);

N = 40;

for i = 1:K
    f{i} = c{i}(1:N, 1:2);
end

%% tomo M de cada clase y calculo las medias iniciales

M = 1;

u = zeros(K, 2);
for i = 1:K
    u(i, :) = mean(f{i}(1:M, :), 1);
end

%% L iteraciones

k = zeros(1, K * N); % supongo que hay N muestras para cada clase

x = zeros(K * N, 2);

% pongo todas las muestras juntas

for i = 1:K
    x((i-1)*N + 1:(i-1)*N + N, :) = f{i}; 
end

%%
L = 15;
f_form = figure;
f_d = figure;
d_iter = zeros(1, L); % distorsion total para cada iteracion

for j = 1:L
    d = zeros(1,K); % distorsion para cada clase
    for i = 1:length(x)
        dist = zeros(1,K);
        for z=1:K
            dist(z) = (x(i, :) - u(z,:)) * (x(i, :) - u(z,:))'; % distancia a la media
        end
        [m,l] = min(dist);
        k(i) = l;
        d(l) = d(l) + m;
    end
    d_iter(j) = sum(d);
 
    xk = cell(1,K);
    
    fprintf('Interacion %i\n',j);
    
    for z=1:K
        xk{z} = x(k==z,:); % x clasificados con la clase z en xk{z}
        fprintf('Total clase %i: %i\n',z,sum(k==z));
    end
    
    fprintf('\n');
    
    % calculo las nuevas medias
    
    for z = 1:K
        u(z, :) = mean(xk{z}, 1);
    end
    
    % graficos
    
    figure(f_form);
    clf;
 
    for z=1:K
        plot(u(z,1),u(z,2),'o','color',col(z),'markersize',10,'markerfacecolor',col(z));
        hold on;
        plot(xk{z}(:, 1), xk{z}(:, 2), 'o','color',col(z));
        hold on;
    end

    xlabel('F1 [Hz]');
    ylabel('F2 [Hz]');
    
    str = sprintf('Iteracion %i', j);
    title(str);
 
    
    figure(f_d);
    clf;
    plot(1:L, d_iter, 'r');
    title(str);
    ylabel('Distorsion');
 
    pause();
end


