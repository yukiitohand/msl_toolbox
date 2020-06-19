function proj_mastcam2MSLDEM_v5_mexw(mastcamdata_obj,MSLDEMdata,MSTprj)
% proj_mastcam2MSLDEM_v5_mexw(mastcamdata_obj,MSLDEMdata,MSTprj)
%   Project mastcam image onto MSLDEMdata
%  INPUTS:
%    mastcamdata_obj: MASTCAMdata or MASTCAMgroup_eye class object
%    MSLDEMdata: HSI class obj, 
%      MSLDEMdata at
%      https://astrogeology.usgs.gov/search/map/Mars/MarsScienceLaboratory/
%      Mosaics/MSL_Gale_DEM_Mosaic_10m
%    MSTprj: object of class MASTCAMCameraProjectionMSLDEM
%
%  OUTPUTS:
%   No outputs. Following properties of MSTprj will be filled.
%     mastcam_NEE: [L_im x S_im x 3]
%       pages 1,2,3 are northing, easting, and elevation. 
%     mastcam_msldemc_ref: [L_im x S_im x 3]
%       indicate which triangle in the MSLDEMdata, the pixel is located.
%       pages 1,2,3 are sample, line, and j. j is either 1 or 2, indicating
%       which triangle at the (sample, line).
%     mastcam_range: [L_im x S_im]
%       range for each pixel.
%     mastcam_msldemc_nn: [L_im x S_im x 2]
%       nearest neighbor pixels in DDR image. The first page is the sample
%       indices and the second is the line indexes.
%     msldemc_imFOVmask_nh: [L_demc x S_demc x 1]
%       Boolean, imFOV_mask with hidden points are removed.
%    

%-------------------------------------------------------------------------%
% Get the size of the mastcam image
%-------------------------------------------------------------------------%
switch class(mastcamdata_obj)
    case 'MASTCAMdata'
        L_im = mastcamdata_obj.hdr.lines; S_im = mastcamdata_obj.samples;
    case 'MASTCAMgroup_eye'
        L_im = mastcamdata_obj.L_im; S_im = mastcamdata_obj.S_im;
    otherwise
        error('second input needs to be either MASTCAMdata or MASTCAMgroup_eye class.');
end

%-------------------------------------------------------------------------%
% Get cam_C_geo
%-------------------------------------------------------------------------%
cmmdl = mastcamdata_obj.CAM_MDL;
rover_nav_coord = mastcamdata_obj.ROVER_NAV;
cmmdl_geo = transform_CAHVOR_MODEL_wROVER_NAV(cmmdl,rover_nav_coord);
C_geo = cmmdl_geo.C;
A_geo = cmmdl_geo.A;

%-------------------------------------------------------------------------%
% Construct Image grid
%-------------------------------------------------------------------------%
imx_im_1d = 0:(S_im-1);
imy_im_1d = reshape(0:(L_im-1),[],1);
[imx_im,imy_im] = meshgrid(imx_im_1d,imy_im_1d);
im_imxy_vec2d = permute(reshape(cat(2,imx_im,imy_im),[L_im*S_im,2]),[2,1,3]);

PmC = cmmdl_geo.get_p_minus_c_from_xy(im_imxy_vec2d);
PmC = reshape(PmC',[L_im,S_im,3]);


%%
%-------------------------------------------------------------------------%
% Get MSL DEM information
%-------------------------------------------------------------------------%
l1 = MSTprj.msldemc_hdr_imxy.line_offset+1;
lend = MSTprj.msldemc_hdr_imxy.line_offset+MSTprj.msldemc_hdr_imxy.lines;
s1 = MSTprj.msldemc_hdr_imxy.sample_offset+1;
send = MSTprj.msldemc_hdr_imxy.sample_offset+MSTprj.msldemc_hdr_imxy.samples;
dem_northing_crop = MSLDEMdata.hdr.y(l1:lend);
dem_easting_crop  = MSLDEMdata.hdr.x(s1:send);

% dem_imFOVd_mask_crop = ~isnan(MSTprj.msldemc_imxy(:,:,1));

%% Main computation.
tic; [im_north,im_east,im_elev,msldem_refx,msldem_refy,msldem_refs,im_range,...
    im_nnx,im_nny] = ...
proj_mastcam2MSLDEM_v4_mex(...
    MSLDEMdata.imgpath,...0
    MSLDEMdata.hdr,...1
    MSTprj.msldemc_hdr_imxy,...2
    dem_northing_crop,...3
    dem_easting_crop,...4
    MSTprj.msldemc_imxy(:,:,1),...5
    MSTprj.msldemc_imxy(:,:,2),...6
    MSTprj.msldemc_imFOVmask,...7
    S_im,L_im,...8,9
    C_geo,A_geo,...10,11
    PmC(:,:,1),PmC(:,:,2),PmC(:,:,3)); toc; % 12,13,14


%% Post computation task.
% Evaluate the neaerest neighbor from the center of each image pixels.
idx_nn = im_nnx*MSTprj.msldemc_hdr_imxy.lines + (im_nny+1);
idx_nn = idx_nn(:);
idx_nn_1d_nisnan = (idx_nn>=1);
idx_nn_1d_ok = idx_nn(idx_nn_1d_nisnan);
img_mask_nh = false(MSTprj.msldemc_hdr_imxy.lines*MSTprj.msldemc_hdr_imxy.samples,1);
img_mask_nh(idx_nn_1d_ok) = true;
img_mask_nh = reshape(img_mask_nh,[MSTprj.msldemc_hdr_imxy.lines,MSTprj.msldemc_hdr_imxy.samples]);

%% Fill propoerties of MASTCAMCameraProjectionMSLDEM object.
MSTprj.mastcam_NEE = cat(3,im_north,im_east,im_elev);
MSTprj.mastcam_msldemc_ref = cat(3,msldem_refx,msldem_refy,msldem_refs);
MSTprj.mastcam_range = im_range;
MSTprj.mastcam_msldemc_nn = cat(3,im_nnx,im_nny);
MSTprj.msldemc_imFOVmask_nh = img_mask_nh;

end
