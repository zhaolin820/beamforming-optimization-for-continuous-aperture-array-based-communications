function [H, Q] = generate_SPD_channel(para, user_loc)
%Calculate the channel of SPDAs
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


d = para.lambda/2; % antenna spacing
Nx = floor((para.Lx(2) - para.Lx(1))/d);
Ny = floor((para.Ly(2) - para.Ly(1))/d);

%% calculate channel of SPDA
H = zeros(Nx*Ny, para.K);
for i = 1:Nx
    for j = 1:Ny
        
        s_x = para.Lx(1) + (i-1/2)*d;
        s_y = para.Ly(1) + (j-1/2)*d;

        for k = 1:para.K
            rk = user_loc(:,k); % location of user k
            [h] = channel_function(para, rk, s_x, s_y);
            H((i-1)*Nx+j,k) = h * sqrt(para.lambda^2/(4*pi)); % Equation (69)
        end
    end
end


%% calculate the channel correlation matrix of SPDA
Q = zeros(para.K, para.K);
for i = 1:para.K
    for j = 1:para.K

        h_i = H(:,i); % channel of user i
        h_j = H(:,j); % channel of user j

        Q(j,i) = h_i'*h_j;
    end
end

end

