classdef MSLGaleMosaicRadius_v3 < MSLGaleDEMMosaic_v3
    
    properties
        OFFSET
    end
     methods
        function obj = MSLGaleMosaicRadius_v3(basename,dirpath,varargin)
            obj@MSLGaleDEMMosaic_v3(basename,dirpath,varargin{:});
            if isfield(obj.lbl.OBJECT_IMAGE,'OFFSET')
                obj.OFFSET = obj.lbl.OBJECT_IMAGE.OFFSET;
            else
                fprintf( ...
                    ['No offset is defined\n', ...
                     'Make sure this is the radius data' ...
                    ]);
            end
        end
        
        
     end

end