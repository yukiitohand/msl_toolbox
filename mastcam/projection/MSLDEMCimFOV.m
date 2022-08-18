classdef MSLDEMCimFOV < handle
    % MSLDEMCimFOV
    %  Storing the data related to image FOV mask (most raw data) 
    %    The properties are supposed to be ENVIRasterMultBandMSLDEMCProj
    %    object, and cropped in the same way. Their sizes are also assumed
    %    to be equal.
    %  Properties
    %    mask     : ENVIRasterMultBandMSLDEMCProj class object
    %      int8 image [msldemc_lines x msldemc_samples]
    %    maskr13  : ENVIRasterMultBandMSLDEMCProj class object
    %      in8 image  [msldemc_lines x msldemc_samples]
    %    pxlctrnn : ENVIRasterMultBandMSLDEMCProj class object
    %      int8 image [msldemc_lines x msldemc_samples]
    %      mask of the image
    %
    properties
        mask
        maskr13
        pxlctrnn
    end
    
    methods
        function obj = MSLDEMCimFOV()
        end
        
    end
end
