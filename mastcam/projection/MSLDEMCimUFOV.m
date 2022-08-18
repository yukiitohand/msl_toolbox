classdef MSLDEMCimUFOV < handle
    % MSLDEMCimUFOV
    %  Storing the data related to image UFOV mask 
    %    The properties are supposed to be ENVIRasterMultBandMSLDEMCProj
    %    object, and cropped in the same way. Their sizes are also assumed
    %    to be equal.
    %  Properties
    %    mask  : ENVIRasterMultBandMSLDEMCProj class object
    %      int8 image [msldemc_lines x msldemc_samples]
    %    xynn  : ENVIRasterMultBandMSLDEMCProj class object
    %      in16 image  [msldemc_lines x msldemc_samples x 2]
    %      The first layer is nearest x integer coordinate in the mastcam
    %      image plane, the second layer is nearest y integer coordinate.
    properties
        mask
        xynn
    end
    
    methods
        function obj = MSLDEMCimUFOV()
        end
        
    end
end
