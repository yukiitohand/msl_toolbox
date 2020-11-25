classdef MASTCAMdataDRCX < MASTCAMdata
    % MASTCAMdataDRCX class
    %  MASTCAM data
    
    properties
        
    end
    
    methods
        function obj = MASTCAMdataDRCX(basename,dirpath,varargin)
            obj@MASTCAMdata(basename,dirpath,varargin{:});
        end
    end
end