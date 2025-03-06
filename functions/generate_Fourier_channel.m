function [H, Q] = generate_Fourier_channel(para, user_loc)
%Calculate the channel based on the Fourier method
%
%  [H, Q] = generate_Fourier_channel(para, user_loc)
%Inputs:
%   para: structure of the initial parameters
%   user_loc: locations of all users
%Outputs:
%   H: approximated wavenumber-domain channel
%   Q: wavenumber-domain channel correlation
%Date: 06/03/2025
%Author: Zhaolin Wang

Lx = para.Lx(2) - para.Lx(1);
Ly = para.Ly(2) - para.Ly(1);

Nx = ceil(Lx/para.lambda); 
Ny = ceil(Ly/para.lambda);

%% calculate the approximated wavenumber-domain channel using the Fourier method
H = zeros((2*Nx+1)*(2*Ny+1), para.K);
for i = 1:(2*Nx+1)
    for j = 1:(2*Ny+1)
       
        n_x = i-1-Nx;
        n_y = j-1-Ny;      
        for k = 1:para.K
            rk = user_loc(:,k); % location of user k

            M = 40; % precision of Gauss-Legendre quadrature   
            [x, w] = GaussLegendre(M); % calculate the coefficients of Gauss-Legendre quadrature  

            h_approx = 0;
            for index = 1:M

                h = channel_function(para, rk, x(index)*Lx/2, x*Ly/2);
                phi = Fourier_basis_function(para, n_x, n_y, x(index)*Lx/2, x*Ly/2); 

                h_approx = h_approx + w(index)* (w.*h).'*phi * (Lx*Ly/4); % Equation (63)
            end

            H((i-1)*(2*Nx+1)+j,k) = h_approx;
        end
    end
end


%% calculate the wavenumber-domain correlation matrix
Q = zeros(para.K, para.K);
for i = 1:para.K
    for j = 1:para.K

        h_i = H(:,i); % channel for user i
        h_j = H(:,j); % channel for user j

        Q(j,i) = h_i'*h_j;
    end
end

end

