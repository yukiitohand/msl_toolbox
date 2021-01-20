classdef MSLGaleDEMMosaic_v3 < MSLGaleMosaic_v3
    
    properties
        
    end
     methods
        function obj = MSLGaleDEMMosaic_v3(basename,dirpath,varargin)
            obj@MSLGaleMosaic_v3(basename,dirpath,varargin{:});
        end
        
        function [subimg] = get_subimage_wPixelRange(obj,xrange,yrange,...
                varargin)
            [subimg] = get_subimage_wPixelRange@MSLGaleMosaic_v3(...
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
        function [elev,xy] = get_elev_wNE(obj,nrthng,estng)
            x = obj.easting2x(estng);
            y = obj.northing2y(nrthng);
            [elev] = obj.get_elev(x,y);
        end
     end

end