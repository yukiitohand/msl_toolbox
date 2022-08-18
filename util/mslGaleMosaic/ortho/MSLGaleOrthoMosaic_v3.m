classdef MSLGaleOrthoMosaic_v3 < MSLGaleMosaic_v3
    
    properties
        
    end
     methods
        function obj = MSLGaleOrthoMosaic_v3(basename,dirpath,varargin)
            obj@MSLGaleMosaic_v3(basename,dirpath,varargin{:});
        end
        
        function [subimg] = get_subimage_wPixelRange(obj,xrange,yrange,...
                varargin)
            [subimg] = get_subimage_wPixelRange@MSLGaleMosaic_v3(...
                obj,xrange,yrange,'Precision','uint8',varargin{:});
        end
     end
end