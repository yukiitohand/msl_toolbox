classdef MASTCAMdataDRCL < MASTCAMdata
    % MASTCAMdataDRCL class
    %  MASTCAM data
    
    properties
        
    end
    
    methods
        function obj = MASTCAMdataDRCL(basename,dirpath,varargin)
            obj@MASTCAMdata(basename,dirpath,varargin{:});
        end
    end
end