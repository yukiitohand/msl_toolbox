function [dalpha,dbeta,dgamma] = mastcam_get_dangle(cmmdl,xy_source,xy_target)
% [dalpha,dbeta,dgamma] = mastcam_get_dangle(cmmdl,xy_source,xy_target)
%  Get rotation angles from the corresponding source and target image
%  pixels.
%  INPUTS
%   cmmdl: Camara model, instance of CAHVOR_MODEL class
%   xy_source: [N x 2] - list of point expressed in the image xy coordinate
%   xy_target: [N x 2] - list of the point expressed in the image xy
%   coordiante
%  OUTPUTS
%   dalpha: rotation along Z axis
%   dbeta : rotation along Y axis
%   dgamma: rotation along X axis

if isempty(cmmdl.Hdash) || isempty(cmmdl.Vdash) || isempty(cmmdl.hs)...
        || isempty(cmmdl.vs) || isempty(cmmdl.hc) || isempty(cmmdl.vc)
    cmmdl.get_image_plane_unit_vectors();
end

N = size(xy_source,1);

% Get source pixel and target pixel vectors.
% pmc_source = cmmdl.get_p_minus_c_from_xy(xy_source');
% pmc_target = cmmdl.get_p_minus_c_from_xy(xy_target');

% evaluate the angle along z axis (yaw in the image)
xy_source_h = [xy_source(:,1) zeros(N,1)];
xy_target_h = [xy_target(:,1) zeros(N,1)];
pmc_source_h = cmmdl.get_p_minus_c_from_xy(xy_source_h');
pmc_target_h = cmmdl.get_p_minus_c_from_xy(xy_target_h');

dalpha = mean(acosd(sum(pmc_source_h.*pmc_target_h,1)));

sgn_alpha = sum( (pmc_target_h - pmc_source_h).*cmmdl.Hdash', 1);
if all(sgn_alpha>0)
    dalpha = -dalpha;
elseif all(sgn_alpha<0)
    % dalpha = dalpha;
else
    error('sgn is not defined');
end

% evaluate the angle along y axis (pitch in the image)
xy_source_v = [zeros(N,1) xy_source(:,2)];
xy_target_v = [zeros(N,1) xy_target(:,2)];
pmc_source_v = cmmdl.get_p_minus_c_from_xy(xy_source_v');
pmc_target_v = cmmdl.get_p_minus_c_from_xy(xy_target_v');

dbeta = mean(acosd(sum(pmc_source_v.*pmc_target_v,1)));

sgn_beta = sum( (pmc_target_v - pmc_source_v).*cmmdl.Vdash', 1);
if all(sgn_beta>0)
    % dbeta = dbeta;
elseif all(sgn_beta<0)
    dbeta = -dbeta;
else
    error('sgn is not defined');
end


% evaluate if rotation component exists or not. not implemented yet.
dgamma = 0;


end