classdef MSTBandpass < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        wavelength
        transmission
    end
    
    methods
        function obj = MSTBandpass(varargin)
            if ~isempty(varargin)
                bandpass = mastcamBandpassRead(varargin{:});
                obj.wavelength   = bandpass.wavelength;
                obj.transmission = bandpass.transmission;
            end
        end
        
        function [bvq] = get_value(obj,wvq,varargin)
            [bvq] = mastcam_getBandpassValue(wvq,obj,varargin{:});
        end
        
        function [mstbp_up] = upsample(obj,varargin)
            mstbp_up = MSTBandpass();
            [bandpass_up] = mastcamBandpassUpsample(obj,varargin{:});
            mstbp_up.wavelength   = bandpass_up.wavelength;
            mstbp_up.transmission = bandpass_up.transmission;
        end
        
        function [spc_conv] = convolve(obj,wv,spc,varargin)
            [spc_conv] = mastcamBandpassConvolve(wv,spc,obj,varargin{:});
        end
        
        function [img_conv] = convolveCRISM(obj,wa,img,varargin)
            [img_conv] = mastcamBandpassConvolveCRISM(wa,img,obj,varargin{:});
        end
        
    end
end