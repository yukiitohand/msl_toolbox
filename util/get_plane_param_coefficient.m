function [plane_param,is_in_face] = get_plane_param_coefficient(ppv1,ppv2,ppv3,pipv,precision,gpu,is_page)
% [plane_param] = get_plane_param_coefficient(ppv1,ppv2,ppv3,tppv,is_gpu)
% Evaluate plane parameters 
%    ppv1 + plane_param(1)(ppv2-ppv1) + plane_param(2)(ppv3-ppv1)
% of the orthogonal projection of a give point vector: pipv onto the plane
% determined by three vectors (ppv1, ppv2, ppv3).
% INPUTS
%  ppv1, ppv2, ppv3: plane positional vectors [3 x 1]
%  pipv: tested positional vector [3 x 1], [3 x N]
%  is_gpu: boolean, whether or not to use GPU.
%  is_page: whether or not to any inputs has the third dimesion.
% OUTPUTS
%  plane_param: paramete vector for a plane
%  is_in_face: whether or not orthogonal projection of pipv fall within the
%  the triangular or not.


pdv1 = ppv2-ppv1; % plane direction vector
pdv2 = ppv3-ppv1;
M = cat(2,pdv1,pdv2);
pipv_inp = pipv - ppv1;

% below opeartion is faster than 
% plane_param = mmx('backslash',M, pipv_inp);
% plane_param = M \ pipv_inp;
if is_page && gpu
    Mpinv = pagefun(@mtimes,pagefun(@mldivide, pagefun(@mtimes,permute(M,[2,1,3]),M), eye(size(M,2),precision,'gpuArray')),permute(M,[2,1,3]));
    plane_param = pagefun(@mtimes,Mpinv,pipv_inp);
elseif is_page && ~gpu
    Mpinv = mmx('backslash',M,repmat(eye(size(M,1)),[1,1,size(M,3)])); 
    plane_param = mmx('mult',Mpinv, pipv_inp);
elseif ~ispage && gpu
    Mpinv = M \ eye(size(M,1),precision,'gpuArray');
    plane_param = Mpinv*pipv_inp;
elseif ~ispage && ~gpu
    Mpinv = M \ eye(size(M,1),precision);
    plane_param = Mpinv*pipv_inp;
end

is_in_face = and( all(plane_param>=0,1), sum(plane_param,1)<=1 );

% Mtmp = [pdv1 pdv2] ./ ldv';
% M = Mtmp(1:2,:)-Mtmp(2:3,:);
% htmp = (ppv1-lopv') ./ ldv';
% h = [-htmp(1)+htmp(2) -htmp(2)+htmp(3)];
% plane_param = M \ h;

end