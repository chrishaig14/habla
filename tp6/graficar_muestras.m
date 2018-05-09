function [ ] = graficar_muestras( xs, clasif, marker ,legends, colors, f1_min, f1_max, f2_min, f2_max, K)
    x_k = cell(1,K);

    for k=1:K
        x_k{k} = xs(clasif==k,:);
    end

    hold on

    for k=1:K
        plot(x_k{k}(:,1),x_k{k}(:,2),marker,'color',colors{k});
        hold on;
    end
    xlim([f1_min, f1_max])
    ylim([f2_min, f2_max])

    legend(legends);
end

