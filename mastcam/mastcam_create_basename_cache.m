function [basename_cache_com] = mastcam_create_basename_cache(mastcamdata_obj,varargin)
% [basename_cache_com] = mastcam_create_basename_cache(mastcamdata_obj)
%  Output common part of the basename of cache files. The basename is in
%  the form that includes side_id, drive_id, pose_id, remote sensing mast
%  motion counter.
%  INPUTS
%   mastcamdata_obj: MASTCAMdata obj or MASTCAMgroup_eye that shares same
%   Rover Navigation model
%  OUTPUTS
%   basename_cache_com: basename of the cache files, common part.
%  OPTIONAL PARAMETERS
%   'ROVER_NAV_VERSION': version of the rover nav
%     (default) >> mastcamdata_obj.ROVER_NAV.version


rover_nav_vr = '';
cam_code = '';
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case {'ROVER_NAV_VERSION'}
                rover_nav_vr = varargin{i+1};
            case 'CAM_CODE'
                cam_code = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

if isempty(rover_nav_vr)
    rover_nav_vr = mastcamdata_obj.ROVER_NAV.version;
end

if iscell(mastcamdata_obj.PRODUCT_ID)
    productID_repre = mastcamdata_obj.PRODUCT_ID{1};
else
    productID_repre = mastcamdata_obj.PRODUCT_ID;
end
propMASTCAMdata = getProp_basenameMASTCAM(productID_repre);

if isempty(cam_code)
    cam_code = propMASTCAMdata.cam_code;
end
if isnumeric(propMASTCAMdata.sol)
    sol = sprintf('%04d',propMASTCAMdata.sol);
end
if isnumeric(propMASTCAMdata.seq_id)
    seq_id = sprintf('%06d',propMASTCAMdata.seq_id);
end
site_id  = mastcamdata_obj.RMC.SITE;
drive_id = mastcamdata_obj.RMC.DRIVE;
pose_id  = mastcamdata_obj.RMC.POSE;
rsm_mc   = mastcamdata_obj.RMC.RSM;
basename_cache_com = sprintf('%s%s%s_site%03ddrive%04dpose%03drsm%03d_%s',sol,cam_code,seq_id,...
    site_id,drive_id,pose_id,rsm_mc,rover_nav_vr);

end