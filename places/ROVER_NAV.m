classdef ROVER_NAV < handle
    %ROVER_NAV class
    %  Properties are the column of telemetry.csv
    %   
    %   version: name of the 
    
    properties
        FRAME
        SITE       % Site ID
        DRIVE      % Drive ID
        POSE       % Pose ID
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
        version
        MAP
        DEM
        RADIUS
        RADIUS_OFFSET
    end
    
    methods
        function obj = ROVER_NAV(varargin)
            obj.FRAME = -1; obj.SITE = -1; obj.DRIVE = -1; obj.POSE = -1;
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
                        case 'VERSION'
                            obj.version = varargin{i+1};
                        case 'MAP'
                            obj.MAP = varargin{i+1};
                        case 'DEM'
                            obj.DEM = varargin{i+1};
                        case 'RADIUS'
                            obj.RADIUS = varargin{i+1};
                        case 'RADIUS_OFFSET'
                            obj.RADIUS_OFFSET = varargin{i+1};
                        otherwise
                            error('Parameter: %s', varargin{i});
                    end
                end
            end
            obj.get_ROVER_NAV_rot_mat();
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
        
        function update_DEM(obj,MSLDEMdata)
            if ~(isa(MSLDEMdata,'MSL_ORBITAL_DEM') ...
                    || isa(MSLDEMdata,'MSLGaleDEMMosaic_v3'))
                error(['The input MSLDEMdata needs to be an object of' ...
                       ' class MSL_ORBITAL_DEM or MSLGaleDEMMosaic_v3']);
            end
            [elev,x,y] = MSLDEMdata.get_elev_wlatlon(...
                obj.LONGITUDE,obj.PLANETOCENTRIC_LATITUDE);
            obj.DEM_PIXEL_LINE   = y-0.5; % from the maximum latitude
            obj.DEM_PIXEL_SAMPLE = x-0.5; % from the westernmost longitude
            obj.ELEVATION = elev;
            obj.DEM = MSLDEMdata.basename;
            if isa(MSLDEMdata,'MSLGaleMosaicRadius_v3')
                obj.RADIUS = elev;
                obj.RADIUS_OFFSET = MSLDEMdata.OFFSET;
            end
        end
        
        function update_MAP(obj,MSLOrthodata)
            if ~(isa(MSLOrthodata,'MSL_ORBITAL_MAP') ...
                    || isa(MSLOrthodata,'MSLGaleOrthoMosaic_v3'))
                error(['The input MSLDEMdata needs to be an object of' ...
                       ' class MSL_ORBITAL_DEM or MSLGaleDEMMosaic_v3']);
            end
            x = MSLOrthodata.lon2x(obj.LONGITUDE);
            y = MSLOrthodata.lat2y(obj.PLANETOCENTRIC_LATITUDE);
            obj.MAP_PIXEL_LINE   = y-0.5; % from the maximum latitude
            obj.MAP_PIXEL_SAMPLE = x-0.5; % from the westernmost longitude
            obj.MAP = MSLOrthodata.basename;
        end
        
    end
end