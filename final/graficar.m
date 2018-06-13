function [] = graficar( x, hmm )
    plot(x(:,1),x(:,2),'or');
    hold on;
    for i=2:4
        plotgaus(hmm.means{i},hmm.vars{i});
        hold on;
    end
end