function [bvq] = mastcam_getBandpassValue(wvq,mstbpqe_obj,varargin)
% [bvq] = mastcam_getBandpassValue(wvq,mstbpqe_obj,varargin)
%   Get the value for bandpass function for give queried wavelength
%  INPUTS
%   wvq: queried wavelength array.
%   mstbpqe_obj: class object of MSTBandpass or MSTQuantumEfficiency 
%                having field "wavelength"
%   propname_bandpass: Name of the property of "mstbpqe_obj" for which 
%      values are obtained. If not specified, "transmission" will be used.
%  OUTPUTS
%   bvq: bandpass value at the queried wavelength
%  
%  USAGE:
%    >> [bvq] = mastcam_getBandpassValue(wvq,bandpass);
%    >> [bvq] = mastcam_getBandpassValue(wvq,bandpass,propname_bandpass);
%   

method = 'linear';
extrap_opt = {};
propname_bandpass = 'transmission';

if isempty(varargin)
    
elseif length(varargin)==1
    propname_bandpass = varargin{1};
else
    error('Input is not right');
end

bvq = interp1(mstbpqe_obj.wavelength,mstbpqe_obj.(propname_bandpass),wvq,method,extrap_opt{:});

end

