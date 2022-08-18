function [plane_param,is_in_face] = get_plane_param_coefficient_2d(ppv1,ppv2,ppv3,pipv)
% [plane_param] = get_plane_param_coefficient(ppv1,ppv2,ppv3,tppv,is_gpu)
% Evaluate plane parameters 
%    ppv1 + plane_param(1)(ppv2-ppv1) + plane_param(2)(ppv3-ppv1)
% of the orthogonal projection of a give point vector: pipv onto the plane
% determined by three vectors (ppv1, ppv2, ppv3).
% INPUTS
%  ppv1, ppv2, ppv3: plane positional vectors [2 x 1]
%  pipv: tested positional vector [2 x 1], [2 x N]
%  is_gpu: boolean, whether or not to use GPU.
%  is_page: whether or not to any inputs has the third dimesion.
% OUTPUTS
%  plane_param: paramete vector for a plane
%  is_in_face: whether or not orthogonal projection of pipv fall within the
%  the triangular or not.


pdv1 = ppv2-ppv1; % plane direction vector
pdv2 = ppv3-ppv1;
% M = cat(2,pdv1,pdv2);
pipv_inp = pipv - ppv1;

Mpinv = [pdv2(2) -pdv2(1); -pdv1(2) pdv1(1)] ./ (pdv1(1)*pdv2(2)-pdv1(2)*pdv2(1));
% Mpinv = M \ eye(size(M,1),precision);
plane_param = Mpinv*pipv_inp;
is_in_face = and( all(plane_param>=0,1), [1,1]*plane_param<=1 );


end