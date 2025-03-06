clc
clear all
close all

addpath("functions");

para = para_init();
num = 100; % number of Monte Carlo simulations
alpha = 1/para.K*ones(para.K,1); % weight of each user

% generate random user location
user_loc_all = zeros(3, para.K, num);
for i=1:num
    user_loc_x = 10*rand(1,para.K)-5;
    user_loc_y = 10*rand(1,para.K)-5;
    user_loc_z = 15*rand(1,para.K)+15;
    
    user_loc_all(:,:,i) = [user_loc_x; user_loc_y; user_loc_z];
end


A = 0.05:0.05:0.5; % aperture size 0.05 m^2 - 0.5 m^2

R_CAPA_all = zeros(length(A), num);
R_CAPA_ZF_all = zeros(length(A), num);
R_Fourier_all = zeros(length(A), num);
R_Fourier_ZF_all = zeros(length(A), num);
R_SPD_all = zeros(length(A), num);
R_SPD_ZF_all = zeros(length(A), num);

WaitMessage = parfor_wait(length(A)*num, 'Waitbar', true);
for i = 1:length(A) 
    para.A = A(i); % m^2
    
    % aperture region
    para.Lx = [-sqrt(para.A)/2,  sqrt(para.A)/2 ];
    para.Ly = [-sqrt(para.A)/2,  sqrt(para.A)/2 ];

    parfor n = 1:num
        user_loc = user_loc_all(:,:,n);

        [Q] = generate_CAPA_channel_correlation(para, user_loc);
        [Rate_CAPA, Rate_CAPA_ZF] = rate_calculator(para, alpha, Q);
        R_CAPA_all(i, n) = Rate_CAPA; R_CAPA_ZF_all(i, n) = Rate_CAPA_ZF;

        [~, Q] = generate_Fourier_channel(para, user_loc);
        [Rate_Fourier, Rate_Fourier_ZF] = rate_calculator(para, alpha, Q);
        R_Fourier_all(i, n) = Rate_Fourier; R_Fourier_ZF_all(i, n) = Rate_Fourier_ZF;

        [~, Q] = generate_SPD_channel(para, user_loc);
        [Rate_SPD, Rate_SPD_ZF] = rate_calculator(para, alpha, Q);
        R_SPD_all(i, n) = Rate_SPD; R_SPD_ZF_all(i, n) = Rate_SPD_ZF;

        WaitMessage.Send;
    end
end
WaitMessage.Destroy



set(groot,'defaultAxesTickLabelInterpreter','latex');
figure; hold on; box on;

h3 = plot(A, mean(R_Fourier_all,2), '-o', 'LineWidth', 1.5, 'Color', '#77AC30');
h4 = plot(A, mean(R_Fourier_ZF_all,2), ':o', 'LineWidth', 1.5, 'Color', '#77AC30');
h1 = plot(A, mean(R_CAPA_all,2), '-sb', 'LineWidth', 1.5);
h2 = plot(A, mean(R_CAPA_ZF_all,2), ':sb', 'LineWidth', 1.5);
h5 = plot(A, mean(R_SPD_all,2), '-dm', 'LineWidth', 1.5, 'Color', '#DE7D00');
h6 = plot(A, mean(R_SPD_ZF_all,2), ':dm', 'LineWidth', 1.5, 'Color', '#DE7D00');

xlabel('Aperture size $A_{\mathrm{T}}$ (m$^2$)', 'Interpreter', 'latex');
ylabel('Weighted Sum-rate (bit/s/Hz)', 'Interpreter', 'latex');
H=gca; H.LineWidth=1;
legend([h1, h2, h3, h4, h5, h6],'CAPA, Proposed method', 'CAPA-ZF, Proposed method',...
    'CAPA, Fourier method', 'CAPA-ZF, Fourier method',...
    'Conventional MIMO', 'Conventional MIMO, ZF', 'Interpreter', 'latex');
ylim([0,10.5]);




function [Rate, Rate_ZF] = rate_calculator(para, alpha, Q)

    [Rate_1] = algorithm_BCD_CoV(para, alpha, Q, 'ZF');
    [Rate_2] = algorithm_BCD_CoV(para, alpha, Q, 'MF');
    Rate = max([Rate_1, Rate_2]);

    [Rate_ZF] = algorithm_ZF(para, alpha, Q);
end