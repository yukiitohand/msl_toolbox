function [mapper_mastcam2msldemc,mapper_msldemc2mastcam_mat,mapcell_msldemc2mastcam] ...
    = mastcam_create_msldemc_mapper(msldemc_imUFOVxynn_obj,mastcam_nn_msldem_obj)
% [mapper_mastcam2msldemc,mapper_msldemc2mastcam_mat,mapcell_msldemc2mastcam] ...
%     = mastcam_create_msldemc_mapper(msldemc_imUFOVxynn_obj,mastcam_nn_msldem_obj)
%  Create a mastcam <-> msldem two-way map projection. Main calculation
%  INPUTS
%   msldemc_imUFOVxynn_obj : ENVIRasterMultBandMSLDEMCProj
%   mastcam_nn_msldem_obj : ENVIRasterMultBand
%  OUTPUTS
%   mapper_mastcam2msldemc     : cell array, {mst_lines,mst_samples}
%      MASTCAM Image pixels -> msldemc image
%   mapper_msldemc2mastcam_mat : matrix [msldemc_lines,msldemc_samples]
%   mapcell_msldemc2mastcam    : cell array, 1-dimensional, its length is
%      equivalent to the number of non-zero components of mapper_msldemc2mastcam_mat
%      msldemc image pixels -> MASTCAM Image pixels
%   
%  OPTIONAL Parameters
%   None
%
%
%%
validateattributes(msldemc_imUFOVxynn_obj, ...
    {'ENVIRasterMultBandMSLDEMCProj'},{},mfilename,'msldemc_imUFOVxynn_obj');
validateattributes(mastcam_nn_msldem_obj, ...
    {'ENVIRasterMultBand'},{},mfilename,'mastcam_nn_msldem_obj');

%%
if isempty(msldemc_imUFOVxynn_obj.img)
    msldemc_imUFOVxynn_img = msldemc_imUFOVxynn_obj.readimg('precision','raw');
else
    msldemc_imUFOVxynn_img = msldemc_imUFOVxynn_obj.img;
end


% pixel correspondence of mastcam_msldemc_nn is updated for the cropped
% image.
if isempty(mastcam_nn_msldem_obj.img)
    mastcam_nn_msldem_img = mastcam_nn_msldem_obj.readimg('precision','raw');
else
    mastcam_nn_msldem_img = mastcam_nn_msldem_obj.img;
end

% converting to the indices to align to msldemc_imUFOVxynn and C-style.
nn_valid = mastcam_nn_msldem_img(:,:,1)>-1;
nn_msldemc_img1 = mastcam_nn_msldem_img(:,:,1);
nn_msldemc_img1(nn_valid) = nn_msldemc_img1(nn_valid) - (msldemc_imUFOVxynn_obj.chdr.sample_offset+1);
nn_msldemc_img2 = mastcam_nn_msldem_img(:,:,2);
nn_msldemc_img2(nn_valid) = nn_msldemc_img2(nn_valid) - (msldemc_imUFOVxynn_obj.chdr.line_offset+1);
mastcam_nn_msldemc_img = cat(3,nn_msldemc_img1,nn_msldemc_img2);

%% Create mapper array
mapper_mastcam2msldemc = msl_create_mapping_mastcam2msldemc_mex_v2(...
    msldemc_imUFOVxynn_img(:,:,1),msldemc_imUFOVxynn_img(:,:,2),...
    mastcam_nn_msldemc_img(:,:,1),mastcam_nn_msldemc_img(:,:,2));

[mapidx_msldemc2mastcam,mapcell_msldemc2mastcam]...
    = msl_create_mapping_msldemc2mastcam_mex_v2(...
        msldemc_imUFOVxynn_img(:,:,1),msldemc_imUFOVxynn_img(:,:,2),...
        mastcam_nn_msldemc_img(:,:,1),mastcam_nn_msldemc_img(:,:,2));

mapper_msldemc2mastcam_mat  = sparse(double(mapidx_msldemc2mastcam+1));

end