function [plane_param_2] = convert_sdtd2st(plane_param,v1,v2,v3,C,A)
% converting the coefficients in the image plane to the coefficient in the
% ground reference coordinate.
% INPUTS
% plane_param: (2 x N)
% v1, v2, v3: (3 x 1), position vector
% C: (3 x 1) position vector of the camera center
% A: (1 x 3) direction vector of the line of sight
% OUTPUTS
% plane_param_2: (2 x N)

dv1 = A*(v1-C);
dv2 = A*(v2-C);
dv3 = A*(v3-C);

pp_num = plane_param ./ [dv2;dv3];
pp_denom = [1 1]*pp_num + (1-[1 1]*plane_param)./dv1;
plane_param_2 = pp_num ./ pp_denom;

end