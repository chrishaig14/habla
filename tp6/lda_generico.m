as = load('a.txt');
os = load('o.txt');
us = load('u.txt');

c = {us,os,as};

legends = ['u';'o';'a'];


%%

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

for i = 1:K
    f{i} = c{i}(1:N, 1:2);
    n_f = n_f + length(f{i});
    t{i} = c{i}(N+1:end,1:2);
    n_t = n_t + length(t{i});
end

%% graficar primeros dos formantes

figure;
for k=1:K
    plot(f{k}(:,1),f{k}(:,2),'o');
    hold on;
end
legend(legends);

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

%% calcular u y sigma

u = zeros(K,2);

for k=1:K
    u(k,:) = mean(f{k},1);
end

sigmas = cell(1,K);

for k=1:K
    sigmas{k} = calcular_sigma(f{k},u(k,:));
end

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
    sigma = sigma + p(k)*sigmas{k};
end


%% test

xs = [];

% pongo todas las muestras juntas

for k = 1:K
    xs = [xs;t{k}];
end

ws = [];

for k=1:K
    ws = [ws , ones(1,length(t{k}))*k];
end

c = zeros(1,length(xs));

for i = 1: length(xs)
    x = xs(i,:);
    % p(k|x)
    p_x = zeros(1,K);
    for k=1:K
        p_x(k) = g_k(x(1),x(2),sigma,u(k,:),p(k));
        
    end
    [m,k_max] = max(p_x);
    
    c(i) = k_max;
end

x_k = cell(1,K);

for k=1:K
    x_k{k} = xs(c==k,:);
end

figure;

for k=1:K
    plot(x_k{k}(:,1),x_k{k}(:,2),'o');
    hold on;
end

legend(legends);


%% calcular error
fprintf('Error: %0.2f %% \n', sum(ws ~= c)/30*100);

f1 = 0:5:2000;
f2 = 0:5:2000;
[F1,F2] = meshgrid(f1,f2);

z = cell(1,K);

fprintf("sigma es: \n");
sigma
u
p


for k=1:K
    z{k} = zeros(size(F1));
    for i=1:numel(F1)
        x = F1(i);
        y = F2(i);
        z{k}(i) = g_k(x,y,sigma,u(k,:), p(k));
    end
    
    figure;
    %size(z{k})
    surf(F1,F2,z{k},'EdgeColor','none');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    view(2)
    str = sprintf('g_k para k=%i',k);
    title(str)
    
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
xlabel('x');
ylabel('y');
zlabel('z');
view(2)
