function [Rate, power_allocation] = algorithm_ZF(para, alpha, Q)
%The proposed corr-ZF algorithm for CAPA beamforming
%
%  [Rate, power_allocation] = algorithm_ZF(para, alpha, Q)
%Inputs:
%   para: structure of the initial parameters
%   alpha: weights of all users
%   Q: channel correlation matrix
%Outputs:
%   R: optimized weighted sum rate
%   power_allocation: optimized power allocation matrix
%Date: 06/03/2025
%Author: Zhaolin Wang

U = inv(Q);
u = real(diag(U)*para.noise); % Equation (52)


%% The following is the classical water filling algorithm. 
% The code refers to https://blog.csdn.net/stay_alive_13/article/details/112846063
% loss: a Kx1 vector
% P: total power to be allocated
% A: diagonal matrix, i.e. weight
% power_allocation: diagonal matrix, power allocation result, trace equals P

loss = u;
A = diag(alpha);
P = para.Pt;

K = length(loss); % number of users

w = diag(A); % width
h = loss./w; % height

allo_set = 1:K; % Initialize the index collection of users to be filled with water
level = (P+sum(loss))/sum(w); % virtual water level
[h_hat, k_hat] = max(h);

while h_hat>=level
    
    allo_set(k_hat) = -1;
    level = (P+sum(loss(allo_set>0)))/sum(w(allo_set>0)); % virtual water level
    [h_hat, ~] = max(h(allo_set>0));
    k_hat = find(h == h_hat);

end

% calculate power allocation matrix
power_allocation = zeros(K,1);
for k = 1:K
    if allo_set(k)>0
        power_allocation(k) = (level - h(k))*w(k);
    end
end


%% calculate weighted sum rate
Rate = sum(alpha.*log2(1 + power_allocation./u)); % Equation (53a)

end

