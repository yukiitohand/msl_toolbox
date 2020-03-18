classdef ROVER_NAV < handle
    %ROVER_NAV class
    %  Properties are the column of telemetry.csv
    %   
    
    properties
        FRAME = -1
        SITE  = -1 % Site ID
        DRIVE = -1 % Drive ID
        POSE  = -1 % Pose ID
        LANDING_X
        LANDING_Y
        LANDING_Z
        NORTHING
        EASTING
        PLANETOCENTRIC_LATITUDE
        PLANETODETIC_LATITUDE
        LONGITUDE
        ELEVATION
        MAP_PIXEL_LINE
        MAP_PIXEL_SAMPLE
        DEM_PIXEL_LINE
        DEM_PIXEL_SAMPLE
        ROLL
        PITCH
        YAW
        SCLK
        SOL
        rot_mat
        rot_mat_inv
    end
    
    methods
        function obj = ROVER_NAV(varargin)
            if (rem(length(varargin),2)==1)
                error('Optional parameters should always go by pairs');
            else
                for i=1:2:(length(varargin)-1)
                    switch upper(varargin{i})
                        case 'FRAME'
                            obj.FRAME = varargin{i+1};
                        case 'SITE'
                            obj.SITE = varargin{i+1};
                        case 'DRIVE'
                            obj.DRIVE = varargin{i+1};
                        case 'POSE'
                            obj.POSE = varargin{i+1};
                        case 'LANDING_X'
                            obj.LANDING_X = varargin{i+1};
                        case 'LANDING_Y'
                            obj.LANDING_Y = varargin{i+1};
                        case 'LANDING_Z'
                            obj.LANDING_Z = varargin{i+1};
                        case 'NORTHING'
                            obj.NORTHING = varargin{i+1};
                        case 'EASTING'
                            obj.EASTING = varargin{i+1};
                        case 'PLANETOCENTRIC_LATITUDE'
                            obj.PLANETOCENTRIC_LATITUDE = varargin{i+1};
                        case 'PLANETODETIC_LATITUDE'
                            obj.PLANETODETIC_LATITUDE = varargin{i+1};
                        case 'LONGITUDE'
                            obj.LONGITUDE = varargin{i+1};
                        case 'ELEVATION'
                            obj.ELEVATION = varargin{i+1};
                        case 'MAP_PIXEL_LINE'
                            obj.MAP_PIXEL_LINE = varargin{i+1};
                        case 'MAP_PIXEL_SAMPLE'
                            obj.MAP_PIXEL_SAMPLE = varargin{i+1};
                        case 'DEM_PIXEL_LINE'
                            obj.DEM_PIXEL_LINE = varargin{i+1};
                        case 'DEM_PIXEL_SAMPLE'
                            obj.DEM_PIXEL_SAMPLE = varargin{i+1};
                        case 'ROLL'
                            obj.ROLL = varargin{i+1};
                        case 'PITCH'
                            obj.PITCH = varargin{i+1};
                        case 'YAW'
                            obj.YAW = varargin{i+1};
                        case 'SCLK'
                            obj.SCLK = varargin{i+1};
                        case 'SOL'
                            obj.SOL = varargin{i+1};
                        otherwise
                            error('Parameter: %s', varargin{i});   
                    end
                end
            end
            get_ROVER_NAV_rot_mat(obj);
        end
        function tf = eq(obj1,obj2)
            if (obj1.SITE == obj2.SITE) ...
                && (obj1.DRIVE == obj2.DRIVE) ...
                && (obj1.POSE == obj2.POSE)
                tf = true;
            else
                tf = false;
            end
        end
        function tf = ne(obj1,obj2)
            if (obj1.SITE ~= obj2.SITE) ...
                || (obj1.DRIVE ~= obj2.DRIVE) ...
                || (obj1.POSE ~= obj2.POSE)
                tf = true;
            else
                tf = false;
            end
        end
        function [] = get_ROVER_NAV_rot_mat(obj)
            obj.rot_mat = get_rot_mat(obj.ROLL,obj.PITCH,obj.YAW);
            obj.rot_mat_inv = get_rot_mat_inv(obj.ROLL,obj.PITCH,obj.YAW);
        end
    end
end