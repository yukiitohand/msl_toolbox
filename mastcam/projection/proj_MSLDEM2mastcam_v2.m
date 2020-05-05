function [MSLDEMprj] = proj_MSLDEM2mastcam_v2(MSLDEMdata,mastcamdata_obj,varargin)
% [MSLDEMprj] = proj_MSLDEM2mastcam_v2(MSLDEMdata,mastcamdata_obj,varargin)
%   evaluate FOV of an image on an ortho-georeferenced image using a
%   georeferenced DEM image.
%  INPUTS:
%    MSLDEMdata: HSI class obj, 
%      MSLDEMdata at
%      https://astrogeology.usgs.gov/search/map/Mars/MarsScienceLaboratory/
%      Mosaics/MSL_Gale_DEM_Mosaic_10m
%    mastcamdata_obj: MASTCAMdata or MASTCAMgroup_eye class object
%  OUTPUTS
%   MSLDEMprj: struct having fields:
%    imFOV_mask: boolean image, [L_dem x S_dem x 1]
%      true if in the FOV, false otherwise.
%      FOV should be true if the linear transformation shows the value
%      between -200 and 200+image_edge. The values outside of this are
%      considered to be inaccurate or invalid. In addition, pixels within
%      the distance 50m are considered to be in FOV, just in case.
%    imxy: [[] x []  x 2] 
%      page 1 is the x coordinate in the image frame,
%      page 2 is the y coordinate in the image frame.
%      x and y coordinate values outside of the range between -200 and 
%      200+image_edge are replaced with nans because the values outside of
%      this is far away from the calibration range, therefore accuracy is
%      not guaranteed. Size depend on FOV. The rectangle that minimally
%      encloses the FOV.

% PROC_MODE = 'naive';
% if (rem(length(varargin),2)==1)
%     error('Optional parameters should always go by pairs');
% else
%     for i=1:2:(length(varargin)-1)
%         switch upper(varargin{i})
%             case 'PROC_MODE'
%                 PROC_MODE = varargin{i+1};
%             otherwise
%                 error('Unrecognized option: %s',varargin{i});
%         end
%     end
% end

%% GET CAMERA, ROVER_NAV, and Image size information

cmmdl = mastcamdata_obj.CAM_MDL;
rover_nav_coord = mastcamdata_obj.ROVER_NAV;

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

%% First compute dem_imFOV_mask
cmmdl_geo = transform_CAHVOR_MODEL_wROVER_NAV(cmmdl,rover_nav_coord);
[dem_imFOV_mask,dem_imFOV_mask_xyd] = get_imFOV_mask_MSLDEM(cmmdl_geo,...
    [L_im,S_im],MSLDEMdata);

%% Next computing image coordinate for each pixel of DEM image.
% This computation is performed after the image mask is calculated because
% until evaluating the FOV, we do not know how to set valid_samples. We
% want to make this less memory intensive. 
[dem_imxy,hdr_dem_imxy] = get_imxy_MSLDEM(cmmdl_geo,MSLDEMdata,...
    dem_imFOV_mask,dem_imFOV_mask_xyd);

%% Summary

MSLDEMprj = [];
MSLDEMprj.imFOV_mask = dem_imFOV_mask;
MSLDEMprj.imxy = dem_imxy;
MSLDEMprj.hdr_imxy = hdr_dem_imxy;

end