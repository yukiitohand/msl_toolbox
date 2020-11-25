classdef MASTCAMdataDRLX < MASTCAMdata
    % MASTCAMdataDRLX class
    %  MASTCAM data
    
    properties
        
    end
    
    methods
        function obj = MASTCAMdataDRLX(basename,dirpath,varargin)
            obj@MASTCAMdata(basename,dirpath,varargin{:});
        end
        
        function [img_iof] = get_IoF(obj)
            if isempty(obj.img)
                obj.readimg();
            end
            switch obj.FILTER_NUMBER
                case 0
                    img_iof = reshape(obj.RADIANCE_FACTOR,[1,1,3]) .* obj.img;
                otherwise
                    img_iof = obj.RADIANCE_FACTOR .* obj.img;
            end
        end
         
    end
end