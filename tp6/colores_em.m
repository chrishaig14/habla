function [ colors_c, F1, F2 ] = colores_em(f1_min, f1_max, f2_min, f2_max,step, u, sigma, p_k, K, colors)
    f1 = f1_min:step:f1_max;
    f2 = f2_min:step:f2_max;
    [F1,F2] = meshgrid(f1,f2);

    gammas = cell(size(F1));

    colors_c = zeros(size(F1,1),size(F1,2),3);

    for i=1:size(F1,1)

                fprintf("%0.1f %% \n",i/size(F1,1)*100);


        for j=1:size(F1,2)

            x = F1(i,j);
            y = F2(i,j);

            gammas{i,j} = zeros(K);

            for k=1:K
                gammas{i,j}(k) = mvnpdf([x,y],u(k,:),sigma{k})*p_k(k);
            end
            p_x = sum(gammas{i,j});
            gammas{i,j} = gammas{i,j}/p_x;

            colors_c(i,j,:) = [0 0 0];
            for k=1:K
                col =  colors{k}*gammas{i,j}(k);
                colors_c(i,j,1) = colors_c(i,j,1) + col(1);
                colors_c(i,j,2) = colors_c(i,j,2) + col(2);
                colors_c(i,j,3) = colors_c(i,j,3) + col(3);
            end
        end

    end

end

