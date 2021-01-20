function [rover_nav] = get_ROVER_NAV_from_mslplc_view_csv(site_id, drive_id, pose_id, varargin)
% [rover_nav] = get_ROVER_NAV_from_mslplc_view_csv(site_id,drive_id,pose_id,varargin)
% Get rover coordinate from a MSLPLC view CSV file for give site_id, 
% drive_id, pose_id. The MSLPLC vidw CSV file can be either 
%   {telemetry.csv, localized_pos.csv,localized_interp.csv,
%    localized_pos_demv2.csv,localized_interp_demv2.csv}.
% INPUTS:
%  site_id: integer, site id
%  drive_id: integer, drive id
%  pose_id: integer, pose id
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
%    'BASENAME_VIEW_CSV': (default) 'telemetry'
%     {'telemetry','localized_pos','localized_interp',...
%      'localized_pos_demv2','localized_interp_demv2'}

global msl_env_vars
localrootDir = msl_env_vars.local_pds_msl_imaging_rootDir;
pds_msl_imaging_URL = msl_env_vars.pds_msl_imaging_URL;

basename_view = 'telemetry';
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'BASENAME_VIEW_CSV'
                basename_view = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

mastcam_rootpath = joinPath(localrootDir,pds_msl_imaging_URL);

dpath_tempfile = '/Users/yukiitoh/src/matlab/toolbox/msl_toolbox/places';

% get telemetry data
mslplc_view_matfname = sprintf('MSLPLC_%s.mat',basename_view);
mslplc_view_matfpath = joinPath(dpath_tempfile,mslplc_view_matfname);
if ~exist(mslplc_view_matfpath,'file')
    lbl_mslplc_view = pds3lblread(joinPath(mastcam_rootpath,'MSLPLC_1XXX/DATA/LOCALIZATIONS',[basename_view, '.lbl']));
    mslplc_view = mslplc_view_csv_read(joinPath(mastcam_rootpath,'MSLPLC_1XXX/DATA/LOCALIZATIONS',[basename_view '.csv']),lbl_mslplc_view);
    MSLPLC_view_FRAME = cat(1,{mslplc_view.FRAME})';
    MSLPLC_view_FRAME = cellfun(@(x) sprintf('% 5s',x),MSLPLC_view_FRAME,'UniformOutput',0);
    MSLPLC_view_FRAME = cat(1,MSLPLC_view_FRAME{:});
    MSLPLC_view_SITE = cat(1,mslplc_view.SITE);
    MSLPLC_view_DRIVE = cat(1,mslplc_view.DRIVE);
    MSLPLC_view_POSE = cat(1,mslplc_view.POSE);
    MSLPLC_view_LANDING_X = cat(1,mslplc_view.LANDING_X);
    MSLPLC_view_LANDING_Y = cat(1,mslplc_view.LANDING_Y);
    MSLPLC_view_LANDING_Z = cat(1,mslplc_view.LANDING_Z);
    MSLPLC_view_NORTHING = cat(1,mslplc_view.NORTHING);
    MSLPLC_view_EASTING = cat(1,mslplc_view.EASTING);
    MSLPLC_view_PLANETOCENTRIC_LATITUDE = cat(1,mslplc_view.PLANETOCENTRIC_LATITUDE);
    MSLPLC_view_PLANETODETIC_LATITUDE = cat(1,mslplc_view.PLANETODETIC_LATITUDE);
    MSLPLC_view_LONGITUDE = cat(1,mslplc_view.LONGITUDE);
    MSLPLC_view_ELEVATION = cat(1,mslplc_view.ELEVATION);
    MSLPLC_view_MAP_PIXEL_LINE = cat(1,mslplc_view.MAP_PIXEL_LINE);
    MSLPLC_view_MAP_PIXEL_SAMPLE = cat(1,mslplc_view.MAP_PIXEL_SAMPLE);
    MSLPLC_view_DEM_PIXEL_LINE = cat(1,mslplc_view.DEM_PIXEL_LINE);
    MSLPLC_view_DEM_PIXEL_SAMPLE = cat(1,mslplc_view.DEM_PIXEL_SAMPLE);
    MSLPLC_view_ROLL = cat(1,mslplc_view.ROLL);
    MSLPLC_view_PITCH = cat(1,mslplc_view.PITCH);
    MSLPLC_view_YAW = cat(1,mslplc_view.YAW);
    MSLPLC_view_SCLK = cat(1,mslplc_view.SCLK);
    MSLPLC_view_SOL = cat(1,mslplc_view.SOL);
    
    save(mslplc_view_matfpath,'MSLPLC_view_FRAME','MSLPLC_view_SITE',...
        'MSLPLC_view_DRIVE','MSLPLC_view_POSE','MSLPLC_view_LANDING_X',...
        'MSLPLC_view_LANDING_Y','MSLPLC_view_LANDING_Z',...
        'MSLPLC_view_NORTHING','MSLPLC_view_EASTING','MSLPLC_view_PLANETOCENTRIC_LATITUDE',...
        'MSLPLC_view_PLANETODETIC_LATITUDE','MSLPLC_view_LONGITUDE',...
        'MSLPLC_view_ELEVATION','MSLPLC_view_MAP_PIXEL_LINE','MSLPLC_view_MAP_PIXEL_SAMPLE',...
        'MSLPLC_view_DEM_PIXEL_LINE','MSLPLC_view_DEM_PIXEL_SAMPLE',...
        'MSLPLC_view_ROLL','MSLPLC_view_PITCH','MSLPLC_view_YAW',...
        'MSLPLC_view_SCLK','MSLPLC_view_SOL');
    
else
    load(mslplc_view_matfpath,'MSLPLC_view_FRAME','MSLPLC_view_SITE',...
        'MSLPLC_view_DRIVE','MSLPLC_view_POSE','MSLPLC_view_LANDING_X',...
        'MSLPLC_view_LANDING_Y','MSLPLC_view_LANDING_Z',...
        'MSLPLC_view_NORTHING','MSLPLC_view_EASTING','MSLPLC_view_PLANETOCENTRIC_LATITUDE',...
        'MSLPLC_view_PLANETODETIC_LATITUDE','MSLPLC_view_LONGITUDE',...
        'MSLPLC_view_ELEVATION','MSLPLC_view_MAP_PIXEL_LINE','MSLPLC_view_MAP_PIXEL_SAMPLE',...
        'MSLPLC_view_DEM_PIXEL_LINE','MSLPLC_view_DEM_PIXEL_SAMPLE',...
        'MSLPLC_view_ROLL','MSLPLC_view_PITCH','MSLPLC_view_YAW',...
        'MSLPLC_view_SCLK','MSLPLC_view_SOL');
    
end

% search matching (site,drive,pose) from the telemetry
loc_mslplc_view = [];
loc_tar = [];
if site_id ~= -1
    loc_mslplc_view = [loc_mslplc_view MSLPLC_view_SITE]; 
    loc_tar = [loc_tar site_id];
end
if drive_id ~= -1
    loc_mslplc_view = [loc_mslplc_view MSLPLC_view_DRIVE]; 
    loc_tar = [loc_tar drive_id];
end
if pose_id ~= -1
    loc_mslplc_view = [loc_mslplc_view MSLPLC_view_POSE]; 
    loc_tar = [loc_tar pose_id];
end 
it = find(all(loc_mslplc_view == loc_tar, 2));

% if nothing matched, then find those that have the same site_id and
% drive_id (ignoring pose_id)
if isempty(it)
    switch basename_view
        case 'telemetry'
            
        case {'localized_interp','localized_interp_demv2','localized_pos','localized_pos_demv2'}
            loc_mslplc_view = [];
            loc_tar = [];
            if site_id ~= -1
                loc_mslplc_view = [loc_mslplc_view MSLPLC_view_SITE]; 
                loc_tar = [loc_tar site_id];
            end
            if drive_id ~= -1
                loc_mslplc_view = [loc_mslplc_view MSLPLC_view_DRIVE]; 
                loc_tar = [loc_tar drive_id];
            end
            it = find(all(loc_mslplc_view == loc_tar, 2));
        otherwise
            error('basename_view = %s is not defined.',basename_view);
    end
end

MSLPLC_view_FRAME_sel = cellstr(MSLPLC_view_FRAME(it,:));
MSLPLC_view_FRAME_sel = strip(MSLPLC_view_FRAME_sel);

if isempty(it)
    rover_nav = ROVER_NAV.empty(0);
else
    for i=1:length(it)
        iit = it(i);
        rover_nav(i) = ROVER_NAV('FRAME',MSLPLC_view_FRAME_sel{i},...
            'SITE',MSLPLC_view_SITE(iit),...
            'DRIVE',MSLPLC_view_DRIVE(iit),...
            'POSE',MSLPLC_view_POSE(iit),...
            'LANDING_X',MSLPLC_view_LANDING_X(iit),...
            'LANDING_Y',MSLPLC_view_LANDING_Y(iit),...
            'LANDING_Z',MSLPLC_view_LANDING_Z(iit),...
            'NORTHING',MSLPLC_view_NORTHING(iit),...
            'EASTING',MSLPLC_view_EASTING(iit),...
            'PLANETOCENTRIC_LATITUDE',MSLPLC_view_PLANETOCENTRIC_LATITUDE(iit),...
            'PLANETODETIC_LATITUDE',MSLPLC_view_PLANETODETIC_LATITUDE(iit),...
            'LONGITUDE',MSLPLC_view_LONGITUDE(iit),...
            'ELEVATION',MSLPLC_view_ELEVATION(iit),...
            'MAP_PIXEL_LINE',MSLPLC_view_MAP_PIXEL_LINE(iit),...
            'MAP_PIXEL_SAMPLE',MSLPLC_view_MAP_PIXEL_SAMPLE(iit),...
            'DEM_PIXEL_LINE',MSLPLC_view_DEM_PIXEL_LINE(iit),...
            'DEM_PIXEL_SAMPLE',MSLPLC_view_DEM_PIXEL_SAMPLE(iit),...
            'ROLL',MSLPLC_view_ROLL(iit),...
            'PITCH',MSLPLC_view_PITCH(iit),...
            'YAW',MSLPLC_view_YAW(iit),...
            'SCLK',MSLPLC_view_SCLK(iit),...
            'SOL',MSLPLC_view_SOL(iit),...
            'version',basename_view,...
            'MAP','msl_orbital_map',...
            'DEM','msl_orbital_dem');
    end
end



end