function [spc_conv] = mastcamBandpassConvolve(wv,spc,mstbpqe_obj,varargin)
% [spc_conv] = mastcamBandpassConvolve(wv,spc,mstbpqe_obj,varargin)
%   Convolve spectra with the bandpass function
% INPUTS
%   wv: wavelength
%   spc: spectrum
%   mstbpqe_obj: class object of MSTBandpass or MSTQuantumEfficiency 
% OUTPUTS
%   spc_conv: convolved value
% OPTIONAL PARAMETERS
%   Any variable length inputs for mstbpqe_obj.get_value

band_inp = mastcam_getBandpassValue(wv,mstbpqe_obj,varargin{:});
% band_inp = mstbpqe_obj.get_value(wv,varargin{:});
spc_conv = sum(spc .* band_inp);

end
