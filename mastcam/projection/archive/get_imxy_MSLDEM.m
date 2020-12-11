function [dem_imxy,hdr_dem_imxy] = get_imxy_MSLDEM(cmmdl,...
    MSLDEMdata,dem_imFOV_mask,dem_imFOV_mask_xyd)
% [dem_imxy,hdr_dem_imxy] = get_imxy_MSLDEM(cmmdl,...
%     MSLDEMdata,dem_imFOV_mask,dem_imFOV_mask_xyd)
%  Obtain imxy coordinate for MSLDEM given imFOV_mask.
%  INPUTS:
%   cmmdl: obj of CAHVOR_MODEL class
%   MSLDEMdata: MSLDEMdata
%   dem_imFOV_mask: [L_dem x S_dem] boolean, fov mask
%   dem_imFOV_mask_xyd: [L_dem x S_dem] boolean, fov mask with valid xyz
%   and direction.
%  OUTPUTS:
%   dem_imxy: boolean (L_dem' x S_dem'), The first page is x-coordinate and
%   the second is y-coordinate.
%   hdr_dem_imxy: storing some offset from the original dem image.

L_dem = MSLDEMdata.hdr.lines; S_dem = MSLDEMdata.hdr.samples;
cmmdl.get_image_plane_unit_vectors();
dem_northing   = reshape(MSLDEMdata.hdr.y,L_dem,1);
dem_easting    = reshape(MSLDEMdata.hdr.x,S_dem,1);

% This computation is performed after the image mask is calculated because
% until evaluating the FOV, we do not know how to set valid_samples. We
% want to make this less memory intensive. 
valid_lines = find(any(dem_imFOV_mask',1));
line_range = valid_lines(1):valid_lines(end);
len_vl = length(line_range);
valid_samples = find(any(dem_imFOV_mask,1));
sample_range = valid_samples(1):valid_samples(end);
len_vs = length(sample_range);

line_offset = line_range(1)-1;
sample_offset = sample_range(1)-1;
% here We made sure that line_range and sample_range are consequtive.

% placeholder for (x,y) coordinate in the image for each pixel of DEM image
dem_imx = nan(len_vs,len_vl);
dem_imy = nan(len_vs,len_vl);

ndem_imFOV_mask_cmmdl_vls = ~dem_imFOV_mask_xyd(line_range,sample_range)';
% pre-asign easting for each line of the data to be analyzed.
deml_g_vs = zeros(len_vs,1);
deml_g_vs(:,2) = dem_easting(sample_range);

MSLDEMdata.fopen_img();
typeName = 'single'; machine = 'ieee-le'; size_type = 4;
fseek(MSLDEMdata.fid_img,size_type*(S_dem*line_offset),-1);

H = cmmdl.H';
V = cmmdl.V';
A = cmmdl.A';
ts_l = size_type*sample_offset;
ts_r = size_type*(S_dem-sample_offset-len_vs);
for li = 1:len_vl
    l = line_range(li);
    %======================================================================
    % Read elevation of the lth line from MSLDEMdata
    %======================================================================
    fseek(MSLDEMdata.fid_img,ts_l,0);
    deml_elev = fread(MSLDEMdata.fid_img,len_vs,typeName,0,machine);
    fseek(MSLDEMdata.fid_img,ts_r,0);
    % deml_elev = MSLDEMdata.lazyEnviReadl(l,0)';
    deml_elev(deml_elev==MSLDEMdata.hdr.data_ignore_value) = nan;
    
    % converting to xyz
    deml_g_vs(:,1) = dem_northing(l);
    deml_g_vs(:,3) = -deml_elev; %(sample_range);
    
    %==============================================================
    % Converting to image coordinate.
    %==============================================================
    deml_g_vec2d_mC = deml_g_vs - cmmdl.C;
    dst2plane_inprod_vec = deml_g_vec2d_mC*A;
    
    deml_imx_vec2d = (deml_g_vec2d_mC*H) ./ dst2plane_inprod_vec;
    deml_imy_vec2d = (deml_g_vec2d_mC*V) ./ dst2plane_inprod_vec;
    

    deml_imx_vec2d(ndem_imFOV_mask_cmmdl_vls(:,li)) = nan;
    deml_imy_vec2d(ndem_imFOV_mask_cmmdl_vls(:,li)) = nan;
    
    dem_imx(:,li) = deml_imx_vec2d;
    dem_imy(:,li) = deml_imy_vec2d;
    
end

MSLDEMdata.fclose_img();

hdr_dem_imxy = MSLDEMdata.hdr;
hdr_dem_imxy.interleave = 'bil';
hdr_dem_imxy.bands = 2;
hdr_dem_imxy.lines = len_vl;
hdr_dem_imxy.samples = len_vs;
hdr_dem_imxy.band_names = {'imx','imy'};
hdr_dem_imxy.line_offset = line_offset;
hdr_dem_imxy.sample_offset = sample_offset;

dem_imx = dem_imx';
dem_imy = dem_imy';
dem_imxy = cat(3,dem_imx,dem_imy);

end