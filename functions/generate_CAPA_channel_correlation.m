function [Q] = generate_CAPA_channel_correlation(para, user_loc)
%Calculate the channel correlation matrix using Gauss-Legendre quadrature 
%
%  [Q] = channel_correlation(para, user_loc)
%Inputs:
%   para: structure of the initial parameters
%   user_loc: locations of all users
%Outputs:
%   Q: channel correlation matrix
%Date: 06/03/2025
%Author: Zhaolin Wang


Q = zeros(para.K, para.K);
for i = 1:para.K
    for j = 1:para.K
        r_i = user_loc(:,i); % location of user i
        r_j = user_loc(:,j); % location of user j
        
        % calculate channel correlation between user i and user j
        q_ij = GuassLegendre_method(para, r_i, r_j);
        Q(j,i) = q_ij;
    end
end


end


function [q_ij] = GuassLegendre_method(para, r_i, r_j)

    M = 40; % precision of Guass-Legendre quadrature

    Lx = para.Lx(2) - para.Lx(1);
    Ly = para.Ly(2) - para.Ly(1);

    % Calculating the channel correlation q_i,j = h_i^H h_j
    [theta, w] = GaussLegendre(M); % calculate the coefficients of Gauss-Legendre quadrature   
    q_ij = 0;
    for index = 1:M       
        H_i = channel_function(para, r_i, theta(index)*Lx/2, theta*Ly/2);
        H_j = channel_function(para, r_j, theta(index)*Lx/2, theta*Ly/2);
        
        q_ij = q_ij + w(index)* (w.*conj(H_i)).'*H_j* (Lx*Ly/4); % Equation (42)
    end

end


