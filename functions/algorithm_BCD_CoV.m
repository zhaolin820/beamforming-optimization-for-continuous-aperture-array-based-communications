function [Rate, W, lambda, mu] = algorithm_BCD_CoV(para, alpha, Q, init_method)
%The proposed BCD-CoV algorithm for CAPA beamforming
%
%  [rate, W, lambda, mu] = algorithm_BCD_CoV(para, alpha, Q, init_method)
%Inputs:
%   para: structure of the initial parameters
%   alpha: weights of all users
%   Q: channel correlation matrix
%   init_method: 'MF' for initialization using matched filering; 'ZF' for initialization using zero forcing
%Outputs:
%   Rate: optimized weighted sum rate
%   [W, lambda, mu]: optimized variables for calculating source current using Equation (27)
%Date: 06/03/2025
%Author: Zhaolin Wang


%% initialization
if strcmp(init_method, 'MF')
    % initialize source current pattern based on matched filtering
    W_init = Q;
elseif strcmp(init_method, 'ZF')
    % initialize source current pattern based on zero forcing
    [~, P] = algorithm_ZF(para, alpha, Q);
    U = inv(Q);
    u = real(diag(U));
    W_init = diag(sqrt(P./u));
else
    error("Wrong specification of initialization method; please specify either 'MF' or 'ZF'.");
end

% initialize the auxiliary variables
mu = zeros(para.K, 1);
lambda = zeros(para.K, 1);
for k = 1:para.K
    w_k = W_init(k,:);
    w_k_sq = abs(w_k).^2;
    mu(k) = sqrt( 1 + w_k_sq(k)/( sum(w_k_sq)- w_k_sq(k) + para.noise ));
    lambda(k) = mu(k)*w_k(k)/(sum(w_k_sq) + para.noise);
end
Rate_pre = sum(alpha.*log2(mu.^2));



%% BCD-CoV iteration
iter_max = 400;
for i = 1:iter_max

    [A, B] = update_A_B(para, mu, lambda, alpha);
    [W] = update_W(para, A, B, Q);
    [mu, lambda] = update_auxiliary_variable(para, A, B, Q, W);
    
    % calculate the rate in the current iteration
    Rate = sum(alpha.*log2(mu.^2));
    
    % check convergence
    if (Rate - Rate_pre)/Rate_pre <= 1e-6
        break;
    end
    Rate_pre = Rate;
end

end


%% Update matrix A and B
function [A, B] = update_A_B(para, mu, lambda, alpha)

A = alpha.*mu.*conj(lambda); % Equation (21)
B = alpha.*abs(lambda).^2; % Equation (21)
C = alpha.*abs(lambda).^2*para.noise/para.Pt; % Equation (21)

A = diag(A./sum(C)); % Equations (26) and (33)
B = diag(B./sum(C)); % Equations (26) and (34)

end

%% Update matrix W
function [W] = update_W(para, A, B, Q)

W = (eye(para.K) + Q*B)\(Q*A); % Equation (37)

end

%% Update auxiliary variables
function [mu, lambda] = update_auxiliary_variable(para, A, B, Q, W)

rho = real(trace( A'*Q*A - 2*real(A*W'*B'*Q) + W'*B'*Q*B*W )); % Equation (38), the transmit power

mu = zeros(para.K, 1);
lambda = zeros(para.K, 1);
for k = 1:para.K
    w_k = W(k,:);
    w_k_square = abs(w_k).^2;
    mu(k) = sqrt( 1 + w_k_square(k)/( sum(w_k_square)- w_k_square(k) + para.noise/para.Pt*rho )); % Equation (39)
    lambda(k) = mu(k)*w_k(k)/(sum(w_k_square) + para.noise/para.Pt*rho); % Equation (40)
end

end