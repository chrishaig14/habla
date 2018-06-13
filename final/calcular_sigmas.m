function sigma = calcular_sigmas(gamma,x,u)

    [numStates, numPts] = size(gamma);
    
    sigma = cell(1,numStates);
    
    for k=2:numStates-1
        sigma{k} = zeros(2,2);
        for n=1:numPts
            sigma{k} = sigma{k} + gamma(k,n)*(x(n,:)-u{k}').'*(x(n,:)-u{k}');
        end
        sigma{k} = sigma{k}/sum(gamma(k,:));
    end
    
    sigma{1} = [];
    sigma{numStates} = [];
end