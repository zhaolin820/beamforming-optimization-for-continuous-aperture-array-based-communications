function [Phi] = Fourier_basis_function(para, n_x, n_y, s_x, s_y)
%Fourier basis function for a point (s_x, s_y, 0) on CAPA
%
%  [Phi] = Fourier_basis_function(para, n_x, n_y, s_x, s_y)
%Inputs:
%   para: structure of the initial parameters
%   [n_x, n_y]: index of the Fourier basis function
%   [s_x, s_y]: x- and y-coordinates of a point on CAPA
%Outputs:
%   Phi: Fourier basis functions
%Date: 06/03/2025
%Author: Zhaolin Wang

Lx = para.Lx(2) - para.Lx(1);
Ly = para.Ly(2) - para.Ly(1);
At = Lx*Ly;

Phi = 1/sqrt(At) * exp(1i*2*pi * ( n_x/Lx.*(s_x - Lx/2) + n_y/Ly.*(s_y - Ly/2) ) );

end

