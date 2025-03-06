function [R, P] = algorithm_WMMSE(para, alpha, H)
%The closed-form WMMSE algorithms for SPDAs developed in [R1]
% [R1] S. S. Christensen, et al, "Weighted sum-rate maximization using weighted MMSE for MIMO-BC beamforming design,"
%      IEEE Transactions on Wireless Communications, vol. 7, no. 12, pp. 4792-4799, December 2008, doi: 10.1109/T-WC.2008.070851.
%
%  [R, P] = algorithm_WMMSE(para, alpha, H)
%Inputs:
%   para: structure of the initial parameters
%   alpha: weights of all users
%   H: channels of all users
%Outputs:
%   R: optimized weighted sum rate rate
%   P: optimized transmit beamformer
%Date: 06/03/2025
%Author: Zhaolin Wang


% initialization
N = size(H, 1);
P = randn(N, para.K) + 1i*randn(N, para.K);

% WMMSE algorithm
R_pre = 0;
for i = 1:20
    [P] = update_P(para, alpha, H, P);

    % check convergence
    [R] = rate(para, alpha, H, P);
    if abs(R - R_pre)/R <= 1e-4
        break;
    end
    R_pre = R;
end

end

%% closed-form WMMSE updates
function [P] = update_P(para, alpha, H, P)

N = size(H, 1);
E = eye(para.K);
Phi = 0; Upsilon = 0;
for k = 1:para.K
    hk = H(:,k);
    pk = P(:,k); 
    I = norm(hk.'*P)^2 + norm(P, 'fro')^2*para.noise/para.Pt; 
    w_k = alpha(k)*(1 + abs(hk.'*pk)^2 / (I - abs(hk.'*pk)^2));
    v_k = hk.'*pk / I;

    Phi = Phi + w_k*abs(v_k)^2 * ( conj(hk)*hk.' + eye(N)*para.noise/para.Pt );
    Upsilon = Upsilon + w_k*conj(v_k)*E(:,k)*hk.';

end

P = Phi'\Upsilon';

end

%% achievable rate
function [R_sum, R] = rate(para, alpha, H, P)

R = zeros(para.K, 1);
for k = 1:para.K
    hk = H(:,k);
    pk = P(:,k); 
    P_I = P; P_I(:,k) = [];
    Ik = norm(hk.'*P_I)^2 + norm(P, 'fro')^2*para.noise/para.Pt; 
    R(k) = log2( 1 + abs(hk.'*pk)^2/Ik );
end
R_sum = sum(alpha.*R);



end

