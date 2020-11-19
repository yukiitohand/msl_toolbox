function [basename_cache_rovnav] = msl_rover_nav_create_basename_cache(mastcamdata_obj,rover_nav_vr)
% [basename_cache_rovnav] = msl_rover_nav_create_basename_cache(mastcamdata_obj,rover_nav_vr)
%  Output common part of the basename of cache files. The basename is in
%  the form that includes side_id, drive_id, pose_id, remote sensing mast
%  motion counter.
%  INPUTS
%   mastcamdata_obj: MASTCAMdata obj or MASTCAMgroup_eye that shares same
%                    Rover Navigation model
%   rover_nav_vr   : version of the rover_nav
%  OUTPUTS
%   basename_cache_rovnav: basename of rover cache file

[basename_cache_com] = mastcam_create_basename_cache(mastcamdata_obj,'ROVER_NAV_VERSION',rover_nav_vr);
basename_cache_rovnav = [basename_cache_com '_ROVER_NAV'];
end