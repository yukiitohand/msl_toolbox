classdef MSTQuantumEfficiency < handle
    % MSTQuantumEfficiency
    %   MASTCAM Bandpass Quantum Efficiency class
    %    Properties: 
    %       wavelength
    %       qe
    %       responsitivity
    %   USAGE:
    %   >> mstqe_red = MSTQuantumEfficiency('Red');
    
    properties
        wavelength
        qe
        responsitivity
    end
    
    methods
        function obj = MSTQuantumEfficiency(varargin)
            if ~isempty(varargin)
                bandpass = mastcamBandpassRead(varargin{:});
                obj.wavelength   = bandpass.wavelength;
                obj.qe = bandpass.qe;
            end
        end
        
        function [bvq] = get_value(obj,wvq,data_type,varargin)
            [bvq] = mastcam_getBandpassValue(wvq,obj,data_type,varargin{:});
        end
        
        function [mstbp_up] = upsample(obj,data_type,varargin)
            mstbp_up = MSTQuantumEfficiency();
            [bandpass_up] = mastcamBandpassUpsample(obj,data_type,varargin{:});
            mstbp_up.wavelength   = bandpass_up.wavelength;
            mstbp_up.(data_type) = bandpass_up.(data_type);
        end
        
        function [spc_conv] = convolve(obj,wv,spc,data_type,varargin)
            [spc_conv] = mastcamBandpassConvolve(wv,spc,obj,data_type,varargin{:});
        end
        
        function [img_conv] = convolveCRISM(obj,wa,img,data_type,varargin)
            [img_conv] = mastcamBandpassConvolveCRISM(wa,img,obj,data_type,varargin{:});
        end
        
        function [R] = qe2responsitivity(obj)
            % respositivity is in the unit of A/W
            hc = 1.23984193; % this value has the uniti of [W/A * um]
            % wavelength is nm.
            R = obj.qe .* obj.wavelength ./ (hc*1000);
            obj.responsitivity = R;
        end
        
    end
end