function [ clasif ] = corregir_etiquetas( ws, c , K)
    %fprintf('Error: %0.2f %% \n', sum(ws ~= c)/length(xs)*100);

    permutaciones = perms(1:K);
    nperms = size(permutaciones,1);
    errores = zeros(1,nperms);

    for n = 1:nperms

        %clasif = zeros(1,length(xs));

        %for i=1:length(xs)

        %             clasif(i) = permutaciones(n,c(i));
        clasif = permutaciones(n,c);

        %end

        errores(n) = sum(ws ~= clasif);

    end

    [error_min, n_min] = min(errores);

    %errores

    %clasif = zeros(1,length(xs));

    %       for i=1:length(xs)

    %clasif(i) = permutaciones(n_min,c(i));
    clasif = permutaciones(n_min,c);

    %      end

    %   nuevo_error = sum(ws ~= clasif);

    %fprintf('Error: %0.2f %% \n', sum(ws ~= clasif)/length(xs)*100);
end

