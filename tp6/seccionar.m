function [ xk, u, theta ] = seccionar( todos, K )
    %% seccionar en K regiones con un angulo theta random

    u = mean(todos,1);
    
    theta = rand()*2*pi/K; % medido desde la horizontal en u

    thetas = [theta:2*pi/K:theta + 2*pi/K*(K-1)];

    thetas = wrapTo2Pi(thetas);

    xk = cell(1,K);

    for i=1:length(todos)
        x_o = todos(i,:);
        x = x_o - u;

        angulo2pi = wrapTo2Pi(atan2(x(2),x(1))); % angulo entre x y u, desde la horizontal

        for k = 1:K-1
            angulo = angulo2pi - thetas(k);
            fin = thetas(k+1) - thetas(k);
            angulo = wrapTo2Pi(angulo);
            fin = wrapTo2Pi(fin);
            if angulo < fin
                xk{k} = [xk{k};x_o];
            end
        end
        angulo = angulo2pi - thetas(K);
        fin = thetas(1) - thetas(K);
        angulo = wrapTo2Pi(angulo);
        fin = wrapTo2Pi(fin);
        if angulo < fin
            xk{K} = [xk{K};x_o];
        end
    end
end

