function [plane_param,is_in_face] = get_plane_param_coefficient(ppv1,ppv2,ppv3,pipv)
% [plane_param] = get_plane_param_coefficient(ppv1,ppv2,ppv3,tppv)
% Evaluate plane parameters 
%    ppv1 + plane_param(1)(ppv2-ppv1) + plane_param(2)(ppv3-ppv1)
% of the orthogonal projection of a give point vector: pipv.
% INPUTS
%  ppv1, ppv2, ppv3: plane positional vectors [3 x 1]
%  tppv: tested positional vector
% OUTPUTS
%  plane_param: paramete vector for a plane
%  is_in_face: whether or not orthogonal projection of pipv fall within the
%  the triangular or not.


pdv1 = ppv2-ppv1; % plane direction vector
pdv2 = ppv3-ppv1;
M = [pdv1 pdv2]; pipv_inp = pipv - ppv1;
plane_param = M \ pipv_inp;
is_in_face = ( all(plane_param>=0) && sum(plane_param)<=1 );

% Mtmp = [pdv1 pdv2] ./ ldv';
% M = Mtmp(1:2,:)-Mtmp(2:3,:);
% htmp = (ppv1-lopv') ./ ldv';
% h = [-htmp(1)+htmp(2) -htmp(2)+htmp(3)];
% plane_param = M \ h;

end