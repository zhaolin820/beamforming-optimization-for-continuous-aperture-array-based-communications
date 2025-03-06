clc
clear all
close all

addpath("functions");

% initialize system parameters
para = para_init();

% generate random user locations
user_loc_x = 10*rand(1,para.K)-5;
user_loc_y = 10*rand(1,para.K)-5;
user_loc_z = 15*rand(1,para.K)+15;
user_loc = [user_loc_x; user_loc_y; user_loc_z];

% weight of each user
alpha = 1/para.K*ones(para.K,1);



%% BCD-CoV algorithm
[Q] = generate_CAPA_channel_correlation(para, user_loc); % channel correlation matrix
init_method='ZF'; % specify initialization method
[Rate_CAPA] = algorithm_BCD_CoV(para, alpha, Q, init_method);
[Rate_CAPA_ZF] = algorithm_ZF(para, alpha, Q);

%% Fourier-based method
[H, Q] = generate_Fourier_channel(para, user_loc);
init_method='ZF'; % specify initialization method
[Rate_Fourier] = algorithm_BCD_CoV(para, alpha, Q, init_method);
[Rate_Fourier_ZF] = algorithm_ZF(para, alpha, Q);
% [Rate_Fourier, P] = algorithm_WMMSE(para, alpha, H); % one can also use the classical WMMSE method 

%% SPDA
[H, Q] = generate_SPD_channel(para, user_loc);
init_method='ZF'; % specify initialization method
[Rate_SPDA] = algorithm_BCD_CoV(para, alpha, Q, init_method);
[Rate_SPDA_ZF] = algorithm_ZF(para, alpha, Q);
% [Rate_Fourier, P] = algorithm_WMMSE(para, alpha, H); % one can also use the classical WMMSE method







