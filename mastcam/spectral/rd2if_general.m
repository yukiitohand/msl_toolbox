function [IoF] = rd2if_general(RD,SF,r)
% [IoF] = rd2if_general(RD,SF,r)
%   convert radiance to I/F
%  Input Parameters
%   RDn: radiance cube [L x S x B]
%   SF: solar spectral irradiance [1 x S x B] at 1AU
%   r : distance from the Sun (AU)
%  Output Parameters
%   IoF : I/F cube [L x S x B]

[L,S,B] = size(RD);

IoF = pi .* RD ./ SF .* (r.^2);


end