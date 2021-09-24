classdef MSLGaleMosaic_v3 < ENVIRasterSingleLayerEquirectProjRot0
    
    properties
        lblpath;
        lbl;
    end
     methods
        function obj = MSLGaleMosaic_v3(basename,dirpath,varargin)
            
            obj@ENVIRasterSingleLayerEquirectProjRot0(...
                basename,dirpath,varargin{:});
            
            obj.lblpath = joinPath(dirpath,[basename '_pds3.lbl']);
            if exist(obj.lblpath,'file')
                obj.lbl = pds3lblread(obj.lblpath);
                obj.hdr = MSLGaleMosaic_v3_lbl2hdr(obj.lbl);
                obj.get_proj_info();
            else
                obj.lblpath = '';
            end

        end
        
        function get_proj_info(obj)
            [obj.proj_info] = mslGaleMosaic_get_cylindrical_proj_info(...
                obj.lbl);
        end
        
     end
end