classdef MASTCAMdataDRLX < MASTCAMdata
    % MASTCAMdataDRLX class
    %  MASTCAM data
    
    properties
        RADIANCE_FACTOR
        RADIANCE_OFFSET
    end
    
    methods
        function obj = MASTCAMdataDRLX(basename,dirpath,varargin)
            obj@MASTCAMdata(basename,dirpath,varargin{:});
            obj.Linearization = 1;
            obj.get_radiance_factor();
        end
        
        function get_radiance_factor(obj)
            [obj.RADIANCE_FACTOR,obj.RADIANCE_OFFSET] = ...
                mastcam_get_radiance_factor(obj.INSTRUMENT_ID,obj.FILTER_NUMBER,obj.lbl);
        end
        
        function [img] = readimg(obj,varargin)
            dtype = 'IoF';
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
            switch upper(dtype)
                case 'RAW'
                    [img] = obj.readimg@MASTCAMdata();
                case {'IOF','IF'}
                    [img] = obj.get_IoF();
                otherwise
                    error('Undefined datatype %s.',dtype);
            end
            if nargout<1
                obj.img = img;
            end
            
        end
        
        function [img_iof] = get_IoF(obj)
            img_raw = obj.readimg('DATATYPE','RAW');
            switch obj.FILTER_NUMBER
                case 0
                    img_iof = reshape(obj.RADIANCE_FACTOR,[1,1,3]) .* img_raw + reshape(obj.RADIANCE_OFFSET,[1,1,3]);
                otherwise
                    img_iof = obj.RADIANCE_FACTOR .* img_raw + obj.RADIANCE_OFFEST;
            end
            % obj.imgIoF = img_iof;
        end
         
    end
end