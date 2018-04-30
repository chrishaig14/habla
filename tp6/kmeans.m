as = load('a.txt');
os = load('o.txt');
us = load('u.txt');

%% mezclar

as = shuffle(as);
os = shuffle(os);
us = shuffle(us);

%% me quedo con N muestras y sÃ³lo los 2 primeros formantes

N = 40;

f_a = as(1:N,1:2);
f_o = os(1:N,1:2);
f_u = us(1:N,1:2);

%% agarro los primeros M y calculo las medias iniciales

M = 1;

u_a = mean(f_a(1:M,:));
u_o = mean(f_o(1:M,:));
u_u = mean(f_u(1:M,:));

%% L iteraciones
c = zeros(1,N);
L = 15;

x = [f_a;f_o;f_u];
 f_form = figure;
 f_d = figure;
 d_iter = zeros(1,L); % distorsión para cada iteración
for j=1:L
    d = [0,0,0]; %distorsion
    for i=1:length(x)
        dist_a = (x(i,:)-u_a)*(x(i,:)-u_a)' % distancia
        dist_o = (x(i,:)-u_o)*(x(i,:)-u_o)'
        dist_u = (x(i,:)-u_u)*(x(i,:)-u_u)'
        [m,l] = min([dist_a,dist_o,dist_u]);
        c(i) = l;
        d(l) = d(l) + m;
    end
    d_iter(j) = sum(d);
    
    x_as = x(c==1,:);
    x_os = x(c==2,:);
    x_us = x(c==3,:);
    
    u_a = mean(x_as);
    u_o = mean(x_os);
    u_u = mean(x_us);
    figure(f_form);
    clf;
    
    x_as = x(c==1,:);
    x_os = x(c==2,:);
    x_us = x(c==3,:);
    colors = ['r*','g*','b*'];

    plot(u_a(1),u_a(2),'ro','markersize',10,'markerfacecolor','r');
    hold on;
    plot(x_as(:,1),x_as(:,2),'ro');
    hold on;
    
    plot(u_o(1),u_o(2),'go','markersize',10,'markerfacecolor','g');
    hold on;
    plot(x_os(:,1),x_os(:,2),'go');
    hold on;
    
    plot(u_u(1),u_u(2),'bo','markersize',10,'markerfacecolor','b');
    hold on;
    plot(x_us(:,1),x_us(:,2),'bo');
    hold on;
    str = sprintf('iteracion %i',j);
    title(str);
    
    figure(f_d);
    clf;
    plot(1:L,d_iter,'r');
    title(str);
    ylabel('distorsion');
    
    pause();
end

x_as = f_a;
x_os = f_o;
x_us = f_u;
colors = ['r*','g*','b*'];

figure;
plot(x_as(:,1),x_as(:,2),'ro');
hold on;
plot(x_os(:,1),x_os(:,2),'go');
hold on;
plot(x_us(:,1),x_us(:,2),'bo');
hold on;

legend('a','o','u');
    

