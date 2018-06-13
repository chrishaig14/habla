function [] = graficar_gammas(x, gamma)
    exp_gamma = exp(gamma);
    [numPts, dim] = size(x);
    for n=1:numPts
        color = exp_gamma(2,n)*[1,0,0];
        color = color + exp_gamma(3,n)*[0,1,0];
        color = color + exp_gamma(4,n)*[0,0,1];
        
        plot(x(n,1),x(n,2),'o','color',color,'markerfacecolor',color);
        hold on
    end

end

