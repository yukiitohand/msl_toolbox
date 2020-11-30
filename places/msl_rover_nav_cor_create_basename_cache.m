function [basename_cache_rovnav] = msl_rover_nav_cor_create_basename_cache(mastcamdata_obj,rover_nav_vr,varargin)
% [basename_cache_rovnav] = msl_rover_nav_cor_create_basename_cache(mastcamdata_obj,rover_nav_vr)
%  Output common part of the basename of cache files. The basename is in
%  the form that includes side_id, drive_id, pose_id, remote sensing mast
%  motion counter.
%  INPUTS
%   mastcamdata_obj: MASTCAMdata obj or MASTCAMgroup_eye that shares same
%                    Rover Navigation model
%   rover_nav_vr   : version of the rover_nav
%  OUTPUTS
%   basename_cache_rovnav: basename of rover cache file

mstcam_code = '';
lr = [];
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            % case {'ROVER_NAV_VERSION'}
            %     rover_nav_vr = varargin{i+1};
            case 'MSTCAM_CODE'
                mstcam_code = varargin{i+1};
            case 'LINEARIZATION'
                lr = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

[basename_cache_com] = mastcam_create_basename_cache(mastcamdata_obj,...
    'ROVER_NAV_VERSION',rover_nav_vr,'CAM_CODE',mstcam_code,'Linearization',lr);
basename_cache_rovnav = [basename_cache_com '_ROVER_NAV'];
end