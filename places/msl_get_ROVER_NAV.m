function [rover_nav,option] = msl_get_ROVER_NAV(mastcamdata_obj,varargin)
% [rover_nav,option] = msl_get_ROVER_NAV(mastcamdata_obj,varargin)
% Get rover coordinate. This is an upper layer function of 
% get_ROVER_NAV_from_mslplc_view_csv.
% INPUTS:
%  mastcamdata_obj: MASTCAMdata obj or MASTCAMgroup_eye obj
% OUTPUTS
%  rover_nav_coord: ROVER_NAV obj
%    having fields: something like
%                       FRAME: 'ROVER'
%                        SITE: 3
%                       DRIVE: 372
%                        POSE: 8
%                   LANDING_X: -0.4610
%                   LANDING_Y: 31.8850
%                   LANDING_Z: 1.2700
%                    NORTHING: -2.7204e+05
%                     EASTING: 8.1468e+06
%     PLANETOCENTRIC_LATITUDE: -4.5895
%       PLANETODETIC_LATITUDE: -4.6438
%                   LONGITUDE: 137.4422
%                   ELEVATION: -4.5020e+03
%              MAP_PIXEL_LINE: 2.1108e+03
%            MAP_PIXEL_SAMPLE: 2.3237e+04
%              DEM_PIXEL_LINE: 528.6200
%            DEM_PIXEL_SAMPLE: 5.8098e+03
%                        ROLL: -1.3700
%                       PITCH: -1.6500
%                         YAW: 116.5000
%                        SCLK: 399622879
%                         SOL: 24
%
% rover frame is referenced from a Site frame:
%  +X: pointing north
%  +Y: east
%  +Z: nadir, pointing down 
% Note that elevation in the CSV files is up positive
% 
%  OPTIONAL Parameters
%    'VERSION': string. 
%       if "VERSION is either of 
%               'telemetry','localized_pos', 'localized_interp',
%               'localized_pos_demv2','localized_interp_demv2',
%       then "get_ROVER_NAV_from_mslplc_view_csv" is called.
%       otherwise, this function try to find a file from a cache folder.
%      
%     

global msl_env_vars
cachedirpath = msl_env_vars.dirpath_cache;

rover_nav_ver = 'localized_interp';
rover_nav_mstcam_code = '';
rover_nav_lr = [];
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case {'VER','VERSION'}
                rover_nav_ver = varargin{i+1};
            case {'MSTCAM_CODE'}
                rover_nav_mstcam_code = varargin{i+1};
            case {'LINEARIZATION'}
                rover_nav_lr = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end



switch lower(rover_nav_ver)
    case {'telemetry','localized_pos','localized_interp','localized_pos_demv2','localized_interp_demv2'}
        site_id  = mastcamdata_obj.RMC.SITE;
        drive_id = mastcamdata_obj.RMC.DRIVE;
        pose_id  = mastcamdata_obj.RMC.POSE;
        [rover_nav] = get_ROVER_NAV_from_mslplc_view_csv(...
            site_id, drive_id, pose_id,'BASENAME_VIEW_CSV',rover_nav_ver);
        option = [];
    otherwise
        [basename_cache_rovnav] = msl_rover_nav_cor_create_basename_cache(...
            mastcamdata_obj,rover_nav_ver,'MSTCAM_CODE',rover_nav_mstcam_code,...
            'LINEARIZATION',rover_nav_lr);
        cachepath = joinPath(cachedirpath,[basename_cache_rovnav '.mat']);
        if exist(cachepath,'file')
            load(cachepath,'rover_nav','option');
        else
            error('File does not exist %s',cachepath);
        end
        
end

end


