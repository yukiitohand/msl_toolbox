classdef MSTSpectralResponse < handle
    % MSTSpectralResponse
    %   MASTCAM Spectral Response class
    %    Properties: 
    %       wavelength
    %       response
    
    properties
        wavelength
        response
    end
    
    methods
        function obj = MSTSpectralResponse(mstqe_obj,mstbp_obj)
            [specResponse] = mastcam_getSpectralResponse(mstqe_obj,mstbp_obj);
            obj.wavelength = specResponse.wavelength;
            obj.response = specResponse.response;

        end
        
        function [bvq] = get_value(obj,wvq,varargin)
            [bvq] = mastcam_getBandpassValue(wvq,obj,'response',varargin{:});
        end
        
        function [mstbp_up] = upsample(obj,varargin)
            mstbp_up = MSTSpectralResponse();
            [bandpass_up] = mastcamBandpassUpsample(obj,'response',varargin{:});
            mstbp_up.wavelength   = bandpass_up.wavelength;
            mstbp_up.response = bandpass_up.response;
        end
        
        function [spc_conv] = convolve(obj,wv,spc,varargin)
            [spc_conv] = mastcamBandpassConvolve(wv,spc,obj,'response',varargin{:});
        end
        
        function [img_conv] = convolveCRISM(obj,wa,img,varargin)
            [img_conv] = mastcamBandpassConvolveCRISM(wa,img,obj,'response',varargin{:});
        end
        
    end
end