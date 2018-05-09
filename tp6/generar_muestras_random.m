function [ todos ] = generar_muestras_random( N, f1_min, f1_max, f2_min, f2_max )
    todos = zeros(N,2);

    for i=1:length(todos)
        todos(i,1) = rand()*(f1_max-f1_min)+f1_min;
        todos(i,2) = rand()*(f2_max-f2_min)+f2_min;
    end
end

