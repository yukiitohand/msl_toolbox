classdef MASTCAM_MSLDEMC_Mapper < handle
    % MASTCAM_MSLDEMC_Mapper
    %   Resolving two-way map projection between MASTCAM and MSLDEM data
    %  
    
    properties
        mastcam2msldemc
        msldemc2mastcam_mat
        mapcell_msldemc2mastcam
        msldem_basename
        msldem_dirpath
        MSLDEMdata
        msldemc_chdr
        msldemc_proj_info
    end
    
    methods
        function obj = MASTCAM_MSLDEMC_Mapper()
            
        end
        
        %% MASK utilities
        function [x_mst,y_mst] = get_mastcam_index(obj,x_msldemc,y_msldemc)
            Lp = length(x_msldemc); x_mst = []; y_mst = [];
            for pi=1:Lp
                if obj.msldemc2mastcam_mat(y_msldemc(pi),x_msldemc(pi))
                    idxes = obj.mapcell_msldemc2mastcam{ ...
                        obj.msldemc2mastcam_mat(y_msldemc(pi),x_msldemc(pi))};
                    x_mst = [x_mst idxes(1,:)];
                    y_mst = [y_mst idxes(2,:)];
                end
            end
        end
        
        function [x_msldemc,y_msldemc] = get_msldemc_index(obj,x_mst,y_mst)
            Lp = length(x_mst); x_msldemc = []; y_msldemc = [];
            for pi=1:Lp
                if obj.mastcam2msldemc{y_mst(pi),x_mst(pi)}
                    idxes = obj.mastcam2msldemc{y_mst(pi),x_mst(pi)};
                    x_msldemc = [x_msldemc idxes(1,:)];
                    y_msldemc = [y_msldemc idxes(2,:)];
                end
            end
        end
        
        % function [selectMask_mastcam,srange,lrange] = get_mastcam_mask_base(obj,x_mst,y_mst)
        %     [selectMask_mastcam,srange,lrange] = create_selectMask(x_mst,y_mst);
        % end
        
        
        function [selectMask_msldemc,srange,lrange] = get_msldemc_mask_base( ...
                obj,x_msldemc,y_msldemc,varargin)
            % [selectMask_msldemc,srange,lrange] = get_msldemc_mask_base(obj,x_mst,y_mst)
            % [selectMask_msldemc,srange,lrange] = get_msldemc_mask_base(obj,x_mst,y_mst,xy_coord)
            %  xy_coord {'NE','LATLON'}
            [selectMask_msldemc,srange,lrange] = create_selectMask(x_msldemc,y_msldemc);
            
            if ~isempty(varargin)
                xy_coord = varargin{1};
                switch upper(xy_coord)
                    case 'PIXEL'
                    case 'NE'
                        srange = obj.msldemc_easting(double(srange));
                        lrange = obj.msldemc_northing(double(lrange));
                    case 'LATLON'
                        srange = obj.msldemc_longitude(double(srange));
                        lrange = obj.msldemc_latitude(double(lrange));
                    otherwise
                        error('Undefined xy_coord %s',xy_coord);
                end
            end
            
        end
        
        function [selectMask_mastcam,srange,lrange] = get_mastcam_mask_base( ...
                obj,x_mst,y_mst)
            [selectMask_mastcam,srange,lrange] = create_selectMask(x_mst,y_mst);
        end
        
        function [selectMask_msldemc,srange,lrange] = get_msldemc_mask( ...
                obj,x_mst,y_mst,varargin)
            % [selectMask_msldemc,srange,lrange] = get_msldemc_mask(obj,x_mst,y_mst,)
            [x_msldemc,y_msldemc] = obj.get_msldemc_index(x_mst,y_mst);
            [selectMask_msldemc,srange,lrange] = get_msldemc_mask_base( ...
                obj,x_msldemc,y_msldemc,varargin{:});
        end
        
        function [selectMask_mastcam,srange,lrange] = get_mastcam_mask( ...
                obj,x_msldemc,y_msldemc)
            [x_mst,y_mst] = obj.get_mastcam_index(x_msldemc,y_msldemc);
            [selectMask_mastcam,srange,lrange] = create_selectMask(x_mst,y_mst);
        end
        
        %% MSLDEMC utilities
        function [lons] = msldemc_longitude(obj,varargin)
            if isempty(varargin)
                lons = obj.msldemc_proj_info.get_longitude(1:obj.msldemc_chdr.samples);
            elseif length(varargin)==1
                lons = obj.msldemc_proj_info.get_longitude(varargin{1});
            else
                error('Invalid input');
            end
        end
        
        function [lats] = msldemc_latitude(obj,varargin)
            if isempty(varargin)
                lats = obj.msldemc_proj_info.get_latitude(1:obj.msldemc_chdr.lines);
            elseif length(varargin)==1
                lats = obj.msldemc_proj_info.get_latitude(varargin{1});
            else
                error('Invalid input');
            end
        end
        
        function [estng] = msldemc_easting(obj,varargin)
            if isempty(varargin)
                estng = obj.msldemc_proj_info.get_easting(1:obj.msldemc_chdr.samples);
            elseif length(varargin)==1
                estng = obj.msldemc_proj_info.get_easting(varargin{1});
            else
                error('Invalid input');
            end
        end
        
        function [nrthng] = msldemc_northing(obj,varargin)
            if isempty(varargin)
                nrthng = obj.msldemc_proj_info.get_northing(1:obj.msldemc_chdr.lines);
            elseif length(varargin)==1
                nrthng = obj.msldemc_proj_info.get_northing(varargin{1});
            else
                error('Invalid input');
            end
        end
        
        function [xi] = msldemc_easting2x(obj,estng)
            xi = obj.msldemc_proj_info.get_x_wEasting(estng);
        end
        function [xi] = msldemc_lon2x(obj,lon)
            xi = obj.msldemc_proj_info.get_x_wlon(lon);
        end
        function [yi] = msldemc_northing2y(obj,nrthng)
            yi = obj.msldemc_proj_info.get_y_wNorthing(nrthng);
        end
        function [yi] = msldemc_lat2y(obj,lat)
            yi = obj.msldemc_proj_info.get_y_wlat(lat);
        end
        
        % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        % lon/easting ctrranges are the ranges of longitude/easting of 
        % [1 hdr.samples]. (the pixel centers of the most left and right
        % pixels)
        % lat/northing ctrranges are the ranges of latitude/northing of 
        % [1 hdr.lines]. (the pixel centers of the upper and lower most
        % pixels)
        function [lat_range] = msldemc_get_lat_ctrrange(obj)
            lat_range = obj.msldemc_latitude([1 obj.msldemc_chdr.lines]);
        end
        function [lon_range] = msldemc_get_lon_ctrrange(obj)
            lon_range = obj.msldemc_longitude([1 obj.msldemc_chdr.samples]);
        end
        function [estng_range] = msldemc_get_easting_ctrrange(obj)
            estng_range = obj.msldemc_easting([1 obj.msldemc_chdr.samples]);
        end
        function [nrthng_range] = msldemc_get_northing_ctrrange(obj)
            nrthng_range = obj.msldemc_northing([1 obj.msldemc_chdr.lines]);
        end
        function [lat_range,lon_range] = msldemc_get_latlon_ctrrange(obj)
            lat_range = obj.msldemc_get_lat_ctrrange();
            lon_range = obj.msldemc_get_lon_ctrrange();
        end
        function [nrthng_range,estng_range] = msldemc_get_NE_ctrrange(obj)
            [estng_range]  = obj.msldemc_get_easting_ctrrange();
            [nrthng_range] = obj.msldemc_get_northing_ctrrange();
        end
        
    end
end
