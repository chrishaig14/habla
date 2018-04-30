%% funciones

function sigma = calcular_sigma(s,u)
    sigma = zeros(2,2);
    for i = 1:length(s)
        x = s(i,:);
        sigma = sigma + (x-u).'*(x-u);
    end
    sigma = sigma * 1/40;
end

