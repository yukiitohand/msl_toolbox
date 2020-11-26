classdef MASTCAMdataDRCX < MASTCAMdata
    % MASTCAMdataDRCX class
    %  MASTCAM data
    
    properties
        
    end
    
    methods
        function obj = MASTCAMdataDRCX(basename,dirpath,varargin)
            obj@MASTCAMdata(basename,dirpath,varargin{:});
        end
        
        function [img] = readimg(obj,varargin)
            dtype = 'uint8';
            if (rem(length(varargin),2)==1)
                error('Optional parameters should always go by pairs');
            else
                for i=1:2:(length(varargin)-1)
                    switch upper(varargin{i})
                        case 'DATATYPE'
                           dtype  = varargin{i+1};
                        otherwise
                            error('Unrecognized option: %s',varargin{i});
                    end
                end
            end
            img = obj.readimg@MASTCAMdata();
            switch upper(dtype)
                case {'UINT','UINT8'}
                    img = uint8(img);
                case {'DOUBLE'}
                    img = double(img);
                otherwise
                    error('Undefined data type %s.',dtype);
            end
            if nargout<1
                obj.img = img;
            end
            
        end
    end
end