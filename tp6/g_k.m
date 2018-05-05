% function [ g ] = g_k( x,y, sigma_k, u_k, pi_k )
%     x = [x,y];
%     g = -1/2*log(det(sigma_k))-1/2*(x-u_k)*inv(sigma_k)*(x-u_k).'+log(pi_k);
% end

function [ g ] = g_k( x,y, det_sigma_k, inv_sigma_k, u_k, pi_k )
    x = [x,y];
    g = -1/2*log(det_sigma_k)-1/2*(x-u_k)*inv_sigma_k*(x-u_k).'+log(pi_k);
end
