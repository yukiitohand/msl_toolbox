function [bpqe_up] = mastcamBandpassUpsample(mstbpqe_obj,varargin)
% [bpqe_up] = mastcamBandpassUpsample(mstbpqe_obj,varargin)
% INPUTS
%   mstbpqe_obj: class object of MSTBandpass or MSTQuantumEfficiency 
%      having field "wavelength"
%   propname_bandpass: Name of the property of "mstbpqe_obj" for which 
%      values are upsampled. If not specified, "transmission" will be used.
% OUTPUTS
%   bpqe_up: struct having two fields, "wavelength" and (propname_bandpass)


propname_bandpass = 'transmission';

if isempty(varargin)
    
elseif length(varargin)==1
    propname_bandpass = varargin{1};
else
    error('Input is not right');
end

wv_step = 1; % [nm]

wv_new = mstbpqe_obj.wavelength(1):wv_step:mstbpqe_obj.wavelength(end);
bpqe_new = interp1(mstbpqe_obj.wavelength,mstbpqe_obj.(propname_bandpass),...
    wv_new,'pchip');

bpqe_up = [];
bpqe_up.wavelength = wv_new;
bpqe_up.(propname_bandpass) = bpqe_new;

end