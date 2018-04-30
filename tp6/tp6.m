as = load('a.txt');
os = load('o.txt');
us = load('u.txt');

%% mezclar

as = shuffle(as);
os = shuffle(os);
us = shuffle(us);

%% me quedo con 40 muestras y sólo los 2 primeros formantes

f_a = as(1:40,1:2);
f_o = os(1:40,1:2);
f_u = us(1:40,1:2);

%% graficar primeros dos formantes

figure;

plot(f_a(:,1),f_a(:,2),'ro');
hold on;

plot(f_o(:,1),f_o(:,2),'bo');
hold on;

plot(f_u(:,1),f_u(:,2),'go');
hold on;

axis equal;

xlabel("F1 [Hz]");
ylabel("F2 [Hz]");

legend('a','o','u');

%% calcular u y sigma

u_a = mean(f_a);
u_o = mean(f_o);
u_u = mean(f_u);

sigma_a = calcular_sigma(f_a,u_a);
sigma_o = calcular_sigma(f_o,u_o);
sigma_u = calcular_sigma(f_u,u_u);

p_a = 1/3;
p_o = 1/3;
p_u = 1/3;

sigma = p_a*sigma_a + p_o*sigma_o + p_u*sigma_u;

%% clasificador

x_a = as(41:50,1:2);
x_o = os(41:50,1:2);
x_u = us(41:50,1:2);

xs = [x_a;x_o;x_u];
ws = [ones(1,10), ones(1,10)*2, ones(1,10)*3];

c = zeros(1,30);

for i = 1: length(xs)
    x = xs(i,:);
    % p(x|a)
    p_x_a = mvnpdf(x, u_a,sigma);
    % p(x|o)
    p_x_o = mvnpdf(x, u_o,sigma);
    % p(x|u)
    p_x_u = mvnpdf(x, u_u,sigma);
    
    [m,l] = max([p_x_a,p_x_o,p_x_u]);
    
    c(i) = l;
end

%% calcular error
fprintf("Error: %0.2f %% \n", sum(ws ~= c)/30*100);

%% graficar elipses de varianza
N = 50,
t = [0:2*pi/(N-1):2*pi];
y1 = sin(t);
y2 = cos(t);

y = [y1;y2];

e_a = chol(sigma_a);
xa = e_a'*y+u_a.';

e_o = chol(sigma_o);
xo = e_o'*y+u_o.';

e_u = chol(sigma_u);
xu = e_u'*y+u_u.';

figure;
plot(xa(1,:),xa(2,:));
hold on;
plot(xo(1,:),xo(2,:));
hold on;
plot(xu(1,:),xu(2,:));
hold on;
size(y)

%figure;

plot(f_a(:,1),f_a(:,2),'ro');
hold on;

plot(f_o(:,1),f_o(:,2),'bo');
hold on;

plot(f_u(:,1),f_u(:,2),'go');
hold on;

x = 0:1500;
y = 0:1500;
[X,Y] = meshgrid(x,y);

xx = [X(:) Y(:)];
p = mvnpdf(xx,u_u,sigma_u);
p = reshape(p,[length(y),length(x)]);
contour(X,Y,p);

axis equal;

xlabel("F1 [Hz]");
ylabel("F2 [Hz]");

legend('a','o','u');

%% funciones

function sigma = calcular_sigma(s,u)
    sigma = zeros(2,2);
    for i = 1:length(s)
        x = s(i,:);
        sigma = sigma + (x-u).'*(x-u);
    end
    sigma = sigma * 1/40;
end

function v = shuffle(a)
    v = a(randperm(length(a)),:);
end