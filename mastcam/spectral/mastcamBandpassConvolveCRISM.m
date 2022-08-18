function [img_conv] = mastcamBandpassConvolveCRISM(wa,img,mstbpqe_obj,varargin)
% [img_conv] = mastcamBandpassConvolveCRISM(wa,img,mstbpqe_obj,varargin)
%   Convolve spectra with the bandpass function
% INPUTS
%   wa: wavelength frame [B x S]
%   img: CRISM image cube [L x S x B]
%   mstbpqe_obj: class object of MSTBandpass or MSTQuantumEfficiency
% OUTPUTS
%   img_conv: convolved image [L x S x 1]
% OPTIONAL PARAMETERS
%   Any variable length inputs for mstbpqe_obj.get_value

[B,S] = size(wa);
[L,Sim,Bim] = size(img);

% If the number of wavelength sampels of wa and img are different, it
% raises an error.
if B~=Bim
    error('Input dimensions of wa and img are not right. See help.');
end

% 'wa' can be either a vector or a matrix. In case of the matrix, wa needs 
% to be have the same numbrer of samples as that of 'img'.
if S>1 && S~=Sim
    error('Input dimensions of wa and img are not right. See help.');
end

band_inp = mastcam_getBandpassValue(wa,mstbpqe_obj,varargin{:});
% band_inp = mstbpqe_obj.get_value(wa,varargin{:});

band_inp = permute(band_inp,[3,2,1]);

img_conv = nansum(img .* band_inp,3);

end