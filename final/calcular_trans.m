function trans = calcular_trans(gamma,xi)
[numStates, numPts] = size(gamma);
nMinOne = numStates - 1;
trans = zeros(numStates, numStates);

for j=2:nMinOne
    for k=2:nMinOne
        
        lognum = logsum(xi(j,k,2:end));  
        logdenom = logsum(gamma(j,1:end-1));
        
%         debería dar lo mismo con esto
%         xi_todo_junto = reshape(xi(j,2:4,2:end),[numel(xi(j,2:4,2:end)),1]);
%         logdenom = logsum(xi_todo_junto);
        
        trans(j,k) = lognum - logdenom;
        
    end
end
    
    % ninguno va al estado 1, del 1 solo se va al 2
    trans(:,1) = log(1e-100);
    trans(1,:) = log(1e-100);
    trans(1,2) = 0;
    
    % ninguno va al estado final, salvo el anterior, y el mismo estado
    % final con prob. 1
    trans(:,numStates) = log(1e-100);
    trans(numStates,1:nMinOne) = log(1e-100);
    trans(numStates,numStates) = 0;
    
    temp = 1 - exp(trans(nMinOne,nMinOne));
    if temp < 1e-100
        temp = 1e-100;
    end
    trans(nMinOne,numStates) = log(temp);
    
    % normalizar
    
    for i=1:numStates
        trans(i,:) = trans(i,:)-logsum(trans(i,:));
    end
    
    trans = exp(trans);
end