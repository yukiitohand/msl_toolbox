classdef MSL_ORBITAL_DEM < HSI
    
    properties
        lblpath;
        lbl;
        proj_info;
    end
     methods
        function obj = MSL_ORBITAL_DEM(basename,dirpath,varargin)
            global msl_env_vars
            localrootDir = msl_env_vars.local_pds_msl_imaging_rootDir;
            pds_msl_imaging_URL = msl_env_vars.pds_msl_imaging_URL;
            [dirpath_guess] = joinPath(localrootDir,pds_msl_imaging_URL,'MSLPLC_1XXX/DATA/MAPS');
            
            % get dirpath if not specified
            if isempty(dirpath)
                dirpath = dirpath_guess;
            end
            if isempty(basename)
               basename = 'msl_orbital_dem'; 
            end
            
            obj@HSI(basename,dirpath,varargin{:});
            obj.lblpath = joinPath(dirpath,[basename '.lbl']);
            obj.lbl = pds3lblread(obj.lblpath);
            obj.hdr = msl_orbital_dem_lbl2hdr(obj.lbl);

        end
        
        function get_proj_info(obj)
            [obj.proj_info] = msl_orbital_dem_cylindrical_proj_info(obj.lbl);
        end
        
        
     end
    
end
    