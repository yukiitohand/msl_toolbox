function [line_param,is_intersect] = line_plane_intersect(lpv1,lpv2,ppv1,ppv2,ppv3)
% [line_param,is_intersect] = line_plane_intersect(lpv1,lpv2,ppv1,ppv2,ppv3)
% Evaluate a line segment determined by lpv1, lpv2 intersects with the
% plane determined by ppv1,ppv2,ppv3.
% INPUTS
%  lpv1, lpv2: line positional vectors [3 x 1]
%  ppv1, ppv2, ppv3: plane positional vectors [3 x 1]
% OUTPUTS
%  line_param: line parameter in an expression lpv1 + t(lpv2-lpv1)
%  is_intersect: whether or not the line segment intersect with the plane

ldv = lpv2-lpv1;
pnv = cross(ppv2-ppv1,ppv3-ppv1); % plane normal vector [3 x 1]
pc = pnv'*ppv1; % plane constant
line_param = (pc-pnv'*lpv1) ./ (pnv'*ldv);
is_intersect = (line_param>=0 && line_param<=1);

end