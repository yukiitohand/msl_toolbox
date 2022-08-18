classdef MSL_ORBITAL_DEM < ENVIRasterSingleLayerEquirectProjRot0
    % MSL_ORBITAL_DEM class
    %  'MSLPLC_1XXX/DATA/MAPS/msl_orbital_dem.img'
    % Usage
    % >> mslorbdem = MSL_ORBITAL_DEM('','');
    % >> 
    properties
        lblpath;
        lbl;
    end
     methods
        function obj = MSL_ORBITAL_DEM(basename,dirpath,varargin)
            global msl_env_vars
            localrootDir = msl_env_vars.local_pds_msl_imaging_rootDir;
            pds_msl_imaging_URL = msl_env_vars.pds_msl_imaging_URL;
            [dirpath_guess] = joinPath(localrootDir,pds_msl_imaging_URL,...
                'MSLPLC_1XXX/DATA/MAPS');
            
            % get dirpath if not specified
            if isempty(dirpath)
                dirpath = dirpath_guess;
            end
            if isempty(basename)
               basename = 'msl_orbital_dem'; 
            end
            
            obj@ENVIRasterSingleLayerEquirectProjRot0(basename,dirpath,...
                varargin{:});
            obj.lblpath = joinPath(dirpath,[basename '.lbl']);
            obj.lbl = pds3lblread(obj.lblpath);
            obj.hdr = msl_places_orbital_map_lbl2hdr(obj.lbl);
            obj.get_proj_info();
            
        end
        
        function get_proj_info(obj)
            [obj.proj_info] = ...
                msl_places_orbital_map_cylindrical_proj_info(obj.lbl);
        end
        
        function [subimg] = get_subimage_wPixelRange(obj,xrange,yrange,...
                varargin)
            [subimg] = get_subimage_wPixelRange@ENVIRasterSingleLayerEquirectProjRot0(...
                obj,xrange,yrange,'Precision','single',varargin{:});
        end
        
        function [elev] = get_elev(obj,x,y)
            [elev] = get_msl_elevation(x,y,obj);
        end
        function [elev,x,y] = get_elev_wlatlon(obj,lon,lat)
            x = obj.lon2x(lon);
            y = obj.lat2y(lat);
            [elev] = obj.get_elev(x,y);
        end
        function [elev,x,y] = get_elev_wNE(obj,nrthng,estng)
            x = obj.easting2x(estng);
            y = obj.northing2y(nrthng);
            [elev] = obj.get_elev(x,y);
        end
     end
    
end
    