function u  = calcular_mus( gamma, x )
    [numStates, numPts] = size(gamma);
    
    u = cell(1,numStates);
    
    for k=2:numStates-1
        u{k} = gamma(k,:)*x/sum(gamma(k,:));
        u{k} = u{k}';
    end
    
    u{1} = [];
    u{numStates} = [];
    
end

