function [para] = para_init()
%Construct a struct of the system parameters 
%  [para] = para_init()
%Inputs:
%   None
%Outputs:
%   para: a struct
%Date: 06/03/2025
%Author: Zhaolin Wang

para.fc = 2.4e9; % carrier frequency 2.4 GHz
c = 3e8; % speed of light
para.lambda = c/para.fc; % wavelength

para.Pt = 0.1; % overall transmit power (A^2)
para.noise = 5.6e-3; % noise power (V^2/m^2)
para.K = 8; % user number

% aperture area
para.A = 0.25; % m^2

% aperture region
para.Lx = [-sqrt(para.A)/2,  sqrt(para.A)/2 ];
para.Ly = [-sqrt(para.A)/2,  sqrt(para.A)/2 ];

% free-space impedance
para.eta = 120*pi;



end

