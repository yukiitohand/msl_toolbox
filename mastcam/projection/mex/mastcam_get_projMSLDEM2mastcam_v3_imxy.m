function [msldemc_imFOVxy] = mastcam_get_projMSLDEM2mastcam_v3_imxy(...
    MSLDEMdata,msldemc_imFOVmask,msldemc_imFOVhdr,L_im,S_im,cmmdl_geo)
% [msldemc_imFOVxy] = mastcam_get_projMSLDEM2mastcam_v3_imxy(...
%     MSLDEMdata,msldemc_imFOVmask,msldemc_imFOVhdr,L_im,S_im,cmmdl_geo)
% computing image coordinate for each pixel of DEM image.
% This computation is performed after the image mask is calculated because
% until evaluating the FOV, we do not know how to set valid_samples. We
% want to make this less memory intensive. 

% L_im = mastcamdata_obj.L_im; S_im = mastcamdata_obj.S_im;

% cmmdl = mastcamdata_obj.CAM_MDL;
% rover_nav_coord = mastcamdata_obj.ROVER_NAV;
% cmmdl_geo = transform_CAHVOR_MODEL_wROVER_NAV(cmmdl,rover_nav_coord);
% cmmdl_geo.get_image_plane_unit_vectors();

msldemc_northing = msldemc_imFOVhdr.y;
msldemc_easting = msldemc_imFOVhdr.x;

tic; [dem_imx,dem_imy] = cahv_get_imxy_MSLDEM_mex(...
    MSLDEMdata.imgpath,MSLDEMdata.hdr,msldemc_imFOVhdr,...
    msldemc_northing,msldemc_easting,...
    msldemc_imFOVmask,S_im,L_im,cmmdl_geo); toc;

msldemc_imFOVxy = cat(3,dem_imx,dem_imy);

end