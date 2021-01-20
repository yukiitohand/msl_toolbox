function [rover_nav_cor,option] = mastcam_correct_rover_nav(rover_nav,cmmdl_rov,MSLDEMdata,varargin)
% [rover_nav_cor] = mastcam_correct_rover_orientation(rover_nav,cmmdl_rov,varargin)
%   Manually correct rover_nav. You can correct rotation and translation
%   with manual inputs. 
% INPUT Parameters
%  rover_nav: object of ROVER_NAV class.
%  cmmdl_rov: object of CAHVORmodel class. Camera model parameters
%             referenced on the Rover coordinate
%  MSLDEMdata: reference DEM data, necessary for getting elevation at the
%              corrected rover location.
% OUTPUT Parameters
%  rover_nav_cor: object of ROVER_NAV class, with location and orientation
%                 fixed.
%  option       : struct, all the optional parameters are stored whether or
%                 not specified.
% OPTIONAL Parameters
%   'ROLL','PITCH','YAW': roll, pitch, yaw angles in degree.
%    (default) 0 for all.
%   'ROTATION_AXIS': {'CAMERA','ROVER','GEOXYZ'}
%     axis along which additional rotation is performed.
%     (default) 'CAMERA'
%     'CAMERA' - the CAMERA coordinate system, X-Y-Z axes are defined as
%                AHV direction. A refers to camera axis vector, 
%                perpendicular to the camera image plane, toward the scene 
%                from the camera center. H is the horizontal vector in the 
%                camera image plane looking right and V is the camera 
%                vertical direction, looking down.
%     'ROVER'  - the Rover Nav coordinate system. X is the rover front
%                direction, 
%     'GEOXYZ' - the coordinate in the base geographical map. In specific, 
%                X-Y-Z direction corresponds to North, East, negative
%                Elevation, respectively.
%   'TRANSLATION_AXIS': {'CAMERA','ROVER','GEOXYZ'}
%     axis along which translation is performed.
%     (default) 'CAMERA'
%   'TRANSLATION_ORDER': {'Before', 'After'}
%     whether or not translation is based on the corrected coordinate 
%     system or not. Valid for the two rotation axes - 'CAMERA' and 'ROVER'.
%     (default) 'Before'
%     'Before' - translation is performed based on the original coordinate
%                system.
%     'After'  - translation is performed based on the corrected coordinate
%                system.
%   'TRANSLATION_VECTOR': 1 x 3 size vector
%     translation parameter for each direction of the coordinate system
%     specified with 'TRANSLATION_AXIS'
%     (default) [0 0 0]

dRolld = 0; dPitchd = 0; dYawd = 0;
rtaxis  = 'CAMERA'; %{'CAMERA','ROVER','GEOXYZ'}
tlaxis  = 'CAMERA';
tlorder = 'BEFORE';
tlvec   = [0 0 0];
vr_stradd  = 'corv1';
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'ROLL'
                dRolld = varargin{i+1};
            case 'PITCH'
                dPitchd = varargin{i+1};
            case 'YAW'
                dYawd = varargin{i+1};
            case 'ROTATION_AXIS'
                rtaxis = upper(varargin{i+1});
            case 'TRANSLATION_AXIS'
                tlaxis = upper(varargin{i+1});
            case 'TRANSLATION_ORDER'
                tlorder = upper(varargin{i+1});
            case 'TRANSLATION_VECTOR'
                tlvec = varargin{i+1};
            case 'VERSION_ADDITIONAL_STRING'
                vr_stradd = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end



% Get camera rotation matrix on the Rover coordinate.
Rcam = [cmmdl_rov.A;cmmdl_rov.Hdash;cmmdl_rov.Vdash]';
% Get camera model and camera rotation matrix
% on the Geographical xyz (north-east-negative elevation) coordiante.
cmmdl_geo = transform_CAHVOR_MODEL_wROVER_NAV(cmmdl_rov,rover_nav);
Rcam_geo = [cmmdl_geo.A;cmmdl_geo.Hdash;cmmdl_geo.Vdash]';

% correcting camera orientation.
switch rtaxis
    case 'CAMERA'
        dR_cam = get_rot_mat(dRolld,dPitchd,dYawd);
        dR_nav = Rcam * dR_cam * Rcam';
        rot_mat_new = rover_nav.rot_mat * dR_nav;
        [rolld,pitchd,yawd] = get_RollPitchYawd_from_rotation_matrix(rot_mat_new);
        
    case 'ROVER'
        % error('Not implemented yet.');
        dR_nav = get_rot_mat(dRolld,dPitchd,dYawd);
        rot_mat_new = rover_nav.rot_mat * dR_nav;
        [rolld,pitchd,yawd] = get_RollPitchYawd_from_rotation_matrix(rot_mat_new);
    case 'GEOXYZ'
        % error('Not implemented yet.');
        rolld = rover_nav.ROLL + dRolld;
        pitchd = rover_nav.PITCH + dPitchd;
        yawd = rover_nav.YAW + dYawd;
        rot_mat_new = get_rot_mat(rolld,pitchd,yawd);
    otherwise
        error('Undefined AXIS %s.',rtaxis);
end

rover_nav_cor = ROVER_NAV('FRAME',rover_nav.FRAME,...
    'SITE',rover_nav.SITE,...
    'DRIVE',rover_nav.DRIVE,...
    'POSE',rover_nav.POSE,...
    'NORTHING',rover_nav.NORTHING,...
    'EASTING',rover_nav.EASTING,...
    'ELEVATION',rover_nav.ELEVATION,...
    'ROLL',rolld,...
    'PITCH',pitchd,...
    'YAW',yawd,...
    'SCLK',rover_nav.SCLK,...
    'SOL',rover_nav.SOL);

% rover_nav_cor = ROVER_NAV('FRAME',rover_nav.FRAME,...
%     'SITE',rover_nav.SITE,...
%     'DRIVE',rover_nav.DRIVE,...
%     'POSE',rover_nav.POSE,...
%     'LANDING_X',rover_nav.LANDING_X,...
%     'LANDING_Y',rover_nav.LANDING_Y,...
%     'LANDING_Z',rover_nav.LANDING_Z,...
%     'NORTHING',rover_nav.NORTHING,...
%     'EASTING',rover_nav.EASTING,...
%     'PLANETOCENTRIC_LATITUDE',rover_nav.PLANETOCENTRIC_LATITUDE,...
%     'PLANETODETIC_LATITUDE',rover_nav.PLANETODETIC_LATITUDE,...
%     'LONGITUDE',rover_nav.LONGITUDE,...
%     'ELEVATION',rover_nav.ELEVATION,...
%     'MAP_PIXEL_LINE',rover_nav.MAP_PIXEL_LINE,...
%     'MAP_PIXEL_SAMPLE',rover_nav.MAP_PIXEL_SAMPLE,...
%     'DEM_PIXEL_LINE',rover_nav.DEM_PIXEL_LINE,...
%     'DEM_PIXEL_SAMPLE',rover_nav.DEM_PIXEL_SAMPLE,...
%     'ROLL',rolld,...
%     'PITCH',pitchd,...
%     'YAW',yawd,...
%     'SCLK',rover_nav.SCLK,...
%     'SOL',rover_nav.SOL);

%% translation
% translation can be performed in 5 different ways.
% You can tranlate the image along three different axes: CAMERA coordinate,
% Rover Nav coordinate, or Geographical XYZ coordinate.
% For the CAMERA and Rover Nav coordinates, you can also choose the current
% system or the corrected system.
switch tlaxis
    case 'CAMERA'
        switch tlorder
            case 'AFTER'
                rot_mat_ref = rot_mat_new * Rcam;
            case 'BEFORE'
                rot_mat_ref = rover_nav.rot_mat * Rcam;
        end
        % rotation
        north_cor = rover_nav_cor.NORTHING + rot_mat_ref(1,:)*tlvec';
        east_cor  = rover_nav_cor.EASTING  + rot_mat_ref(2,:)*tlvec';
    case 'ROVER'
        switch tlorder
            case 'AFTER'
                rot_mat_ref = rot_mat_new;
            case 'BEFORE'
                rot_mat_ref = rover_nav.rot_mat;
        end
        north_cor = rover_nav_cor.NORTHING + rot_mat_ref(1,:)*tlvec';
        east_cor  = rover_nav_cor.EASTING  + rot_mat_ref(2,:)*tlvec';
    case 'GEOXYZ'
        north_cor = rover_nav_cor.NORTHING + tlvec(1);
        east_cor  = rover_nav_cor.EASTING  + tlvec(2);
    otherwise
        error('Undefined AXIS %s.',tlaxis);
end

x = MSLDEMdata.easting2x(east_cor);
y = MSLDEMdata.northing2y(north_cor);
[el_cor]  = get_msl_elevation(x,y,MSLDEMdata);
rover_nav_cor.DEM_PIXEL_LINE = y-0.5;
rover_nav_cor.DEM_PIXEL_SAMPLE = x-0.5;
rover_nav_cor.PLANETOCENTRIC_LATITUDE = MSLDEMdata.latitude(y);
rover_nav_cor.LONGITUDE = MSLDEMdata.longitude(x);
rover_nav_cor.DEM = MSLDEMdata.basename;
rover_nav_cor.NORTHING  = north_cor;
rover_nav_cor.EASTING   = east_cor ;
rover_nav_cor.ELEVATION = el_cor   ;

rover_nav_cor.version = [rover_nav.version '_' vr_stradd];

%%
option = struct('ROTATION_AXIS',rtaxis,'ROLLd',dRolld,'PITCHd',dPitchd,'YAWd',dYawd,...
    'TRANSLATION_AXIS',tlaxis,'TRANSLATION_ORDER',tlorder,'TRANSLATION_VECTOR',tlvec);


end