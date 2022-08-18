function [msldemc_imUFOVmask_crop,msldemc_imUFOVhdr]...
    = mastcam_crop_msldemc_imUFOVmask_v2(msldemc_imFOVmask_obj,msldemc_imUFOVmask)
% [msldemc_imUFOVmask_crop,msldemc_imUFOVhdr]...
%     = mastcam_crop_msldemc_imUFOVmask_v2(msldemc_imFOVmask_obj,msldemc_imUFOVmask)
% Cut the margin of imUFOVmask 
%  INPUTS
%    msldemc_imFOVmask_obj : object of ENVIRasterSingleLayerMSLDEMCProj
%    msldemc_imUFOVmask    : int8 matrix with the size of 
%       [msldemc_imFOVmask.hdr.lines,msldemc_imFOVmask.hdr.samples]
%  OUTPUTS
%    msldemc_imUFOVmask_crop : int8 matrix with the size of 
%       [msldemc_imUFOVhdr.lines, msldemc_imUFOVhdr.samples]
%    msldemc_imUFOVhdr       : struct, having the following four fields
%      samples,lines, sample_offset, line_offset
%      
% 

% 
validateattributes(msldemc_imFOVmask_obj, ...
    {'ENVIRasterSingleLayerMSLDEMCProj','ENVIRasterMultBandMSLDEMCProj'},  ...
    {},mfilename,'msldemc_imFOVmask_obj');
validateattributes(msldemc_imUFOVmask,{'int8'}, ...
    {'size', [msldemc_imFOVmask_obj.hdr.lines,msldemc_imFOVmask_obj.hdr.samples]}, ...
    mfilename,'msldemc_imUFOVmask');

% Crop the mask
valid_lines   = find(any(msldemc_imUFOVmask',1));
lrnge         = [valid_lines(1), valid_lines(end)];
len_vl        = lrnge(2)-lrnge(1)+1;
valid_samples = find(any(msldemc_imUFOVmask,1));
srnge         = [valid_samples(1), valid_samples(end)];
len_vs        = srnge(2)-srnge(1)+1;

% Crop the mask
msldemc_imUFOVmask_crop = msldemc_imUFOVmask(lrnge(1):lrnge(2),srnge(1):srnge(2));

msldemc_imUFOVhdr = [];
msldemc_imUFOVhdr.lines   = len_vl;
msldemc_imUFOVhdr.samples = len_vs;
msldemc_imUFOVhdr.line_offset   = msldemc_imFOVmask_obj.chdr.line_offset+lrnge(1)-1;
msldemc_imUFOVhdr.sample_offset = msldemc_imFOVmask_obj.chdr.sample_offset+srnge(1)-1;

end
