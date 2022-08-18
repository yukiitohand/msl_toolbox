function [mask_hsi,mask_msldemc_UFOVmask_hsiresol,mask_mastcam_hsiresol]...
    = mastcam_get_HSISelectMask_from_msldemcUFOVmask(...
    MSTproj,objHSIdata,x_east,y_north,mask_msldemc_UFOVmask)
% [mask_hsi,mask_msldemc_UFOVmask_hsiresol,mask_mastcam_hsiresol]...
%     = mastcam_get_HSIdataGLTproj_SelectMask_from_msldemcUFOVmask(...
%     MSTproj,objHSIdataGLTproj,x_east,y_north,mask_msldemc_UFOVmask)
% get msldemc_UFOVmask at hsi resolution and 
% INPUT parameters
%  MSTproj: object of MASTCAMCameraProjection
%  objHSIdata: object of HSIdata, needs to be geographically
%  referenced with planetocenteric coordinate system. For HSIdataGLTproj
%  object, pass its GLTdata.
%  x_east: selected easting coordinate point
%  y_north: selected northing coordiante point.
%  mask_msldemc_UFOVmask: associated mask_msldemc_UFOVmask for which a hsi
%  resolution is obtained
% OUTPUT Parameters
%  mask_hsi: mask at the resolution of hsi image
%  msldemc_UFOVmask_hsiresol:
%  mastcam_hsiresol: mask at the resolution of hsi image projected onto
%  MASTCAM image space.

S_hsi = objHSIdata.hdr.samples;
L_hsi = objHSIdata.hdr.lines;

% if isa(objHSIdataGLTproj,'HSIdataGLTproj')
%     S_hsi = objHSIdataGLTproj.GLTdata.hdr.samples;
%     L_hsi = objHSIdataGLTproj.GLTdata.hdr.lines;
% elseif isa(objHSIdataGLTproj,'HSI')
%     S_hsi = objHSIdataGLTproj.hdr.samples;
%     L_hsi = objHSIdataGLTproj.hdr.lines;
% end

mask_hsi = false(L_hsi,S_hsi);
% nearest HSI image pixel.
[xi,yi] = objHSIdata.get_xy_fromNE(x_east,y_north);
% if isa(objHSIdataGLTproj,'HSIdataGLTproj')
%     [xi,yi] = objHSIdataGLTproj.get_GLTxy_fromNE(x_east,y_north);
% elseif isa(objHSIdataGLTproj,'HSI')
%     [xi,yi] = objHSIdataGLTproj.get_xy_fromNE(x_east,y_north);
% end
mask_hsi(yi,xi) = true;


[pos_row,pos_col] = find(mask_msldemc_UFOVmask);
for i=1:length(pos_row)
    rowi = pos_row(i); coli = pos_col(i);
    [xei,yni] = MSTproj.get_NE_from_msldemc_imUFOVxy(coli,rowi);
    % Get the HSI pixel nearest from one selected pixel in
    % msldemc_UFOVmask.
    [xi,yi] = objHSIdata.get_xy_fromNE(xei,yni);
%     if isa(objHSIdataGLTproj,'HSIdataGLTproj')
%         [xi,yi] = objHSIdataGLTproj.get_GLTxy_fromNE(x_east,y_north);
%     elseif isa(objHSIdataGLTproj,'HSI')
%         [xi,yi] = objHSIdataGLTproj.get_xy_fromNE(x_east,y_north);
%     end
    % mask_hsi_ar(i,1) = xi; mask_hsi_ar(i,2) = yi; 
    mask_hsi(yi,xi) = true;
end


mask_msldemc_UFOVmask_hsiresol = false(MSTproj.msldemc_imUFOVhdr.lines,...
    MSTproj.msldemc_imUFOVhdr.samples);
mask_mastcam_hsiresol = false(MSTproj.MASTCAMdata.L_im,MSTproj.MASTCAMdata.S_im);
[hsirow,hsicol] = find(mask_hsi);
for i=1:length(hsirow) 
    hsirowi = hsirow(i); hsicoli = hsicol(i);
    [xrnge_east,yrnge_north] = objHSIdata.get_pixel_rangeNE_fromxy(hsicoli,hsirowi);
%     if isa(objHSIdataGLTproj,'HSIdataGLTproj')
%         [xrnge_east,yrnge_north] = objHSIdataGLTproj.get_pixel_rangeNE_fromGLTxy(hsicoli,hsirowi);
%     elseif isa(objHSIdataGLTproj,'HSI')
%         [xrnge_east,yrnge_north] = objHSIdataGLTproj.get_pixel_rangeNE_fromxy(hsicoli,hsirowi);
%     end
    [msldemcUFOV_hsii_mask,mask_mastcam_i] = MSTproj.get_rangeUFOVmask(xrnge_east,yrnge_north);
    mask_msldemc_UFOVmask_hsiresol = or(mask_msldemc_UFOVmask_hsiresol,msldemcUFOV_hsii_mask);
    mask_mastcam_hsiresol = or(mask_mastcam_hsiresol,mask_mastcam_i);
end
end