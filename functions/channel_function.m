function [H] = channel_function(para, r, s_x, s_y)
%Calculate the channel between a point (s_x, s_y, 0) on CAPA and the user located at r
%
%  [H] = channel_function(para, r, s_x, s_y)
%Inputs:
%   para: structure of the initial parameters
%   r: location vector of the user
%   [s_x, s_y]: x- and y-coordinates of a point on CAPA
%Outputs:
%   H: channel
%Date: 06/03/2025
%Author: Zhaolin Wang


z = sqrt( (r(1,:) - s_x).^2 + (r(2,:) - s_y).^2 + r(3,:).^2 );
H = - 1i*para.eta*exp(-1i*2*pi/para.lambda*z)./(2*para.lambda*abs(z))...
    .*(1 -  (r(2,:) - s_y).^2 ./ (z.^2)); % Equation (8)

end

