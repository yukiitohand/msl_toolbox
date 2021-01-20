function [elev_tar] = get_msl_elevation(x,y,MSLDEMdata)
% [elev_tar] = get_msl_elevation(x,y,MSLDEMdata)
%   Get the elevation at image pixel (x,y) by triangulation of the dem
%   pixel points
% INPUTS
%  x: scalar, x-coordinate
%  y: scalar, y-coordinate
%  MSLDEMdata: MSLGaleDEMMosaic_v3 or MSL_ORBITAL_DEM
% OUTPUTS
%  elev_tar: scalar, elevation at (x,y)

% northing easting target vector
ne_tar = [MSLDEMdata.northing(y);MSLDEMdata.easting(x)];

% left upper pixel of the square region.
yflr = floor(y); xflr = floor(x);


% northing, easting, elevation of the square region
nrthng_llp1 = MSLDEMdata.northing([yflr yflr+1]);
estng_ssp1  = MSLDEMdata.easting([xflr xflr+1]);
msldem_elev = MSLDEMdata.get_subimage_wPixelRange(...
    [xflr xflr+1],[yflr yflr+1],'precision','double');

for j=1:2
    if j==1
        % plane position vector in image space
        ppv1 = [nrthng_llp1(1);estng_ssp1(1)];
        ppv2 = [nrthng_llp1(1);estng_ssp1(2)];
        ppv3 = [nrthng_llp1(2);estng_ssp1(1)];
        el1 = msldem_elev(1,1);
        el2 = msldem_elev(1,2);
        el3 = msldem_elev(2,1);
    elseif j==2
        ppv1 = [nrthng_llp1(1);estng_ssp1(2)]; 
        ppv2 = [nrthng_llp1(2);estng_ssp1(2)];
        ppv3 = [nrthng_llp1(2);estng_ssp1(1)];
        el1 = msldem_elev(1,2);
        el2 = msldem_elev(2,2);
        el3 = msldem_elev(2,1);
    end
    [plane_param,is_in_face] = get_plane_param_coefficient_2d(...
        ppv1,ppv2,ppv3,ne_tar);
    
    if is_in_face
        elev_tar= el1 + ([el2 el3] - el1)*plane_param;
    end
end



end