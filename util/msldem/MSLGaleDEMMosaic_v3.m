classdef MSLGaleDEMMosaic_v3 < HSI
    
    properties
        lblpath;
        lbl;
    end
     methods
        function obj = MSLGaleDEMMosaic_v3(basename,dirpath,varargin)
            
            obj@HSI(basename,dirpath,varargin{:});
            obj.lblpath = joinPath(dirpath,[basename '.lbl']);
            obj.lbl = pds3lblread(obj.lblpath);
            obj.hdr = MSLGaleDEMMosaic_v3_lbl2hdr(obj.lbl);

        end
     end
    
end