function [dem_imFOV_mask,dem_imFOV_mask_xyd] = get_imFOV_mask_MSLDEM(cmmdl,camim_size,MSLDEMdata)
% [dem_imFOV_mask,dem_imFOV_mask_xyd] = get_imFOV_mask_MSLDEM(cmmdl,camim_size,MSLDEMdata)
%  Obtain image FOV for MSLDEM given camera model 
% INPUTS:
%   cmmdl: obj of CAHVOR_MODEL class
%   camim_size: size of the camera image [L_im, S_im]
%   MSLDEMdata: 
%  OUTPUTS:
%   dem_imFOV_mask: boolean (L_dem x S_dem), if a pixel is potentially in FOV,
%   then true.
%   dem_imFOV_mask_xyd: boolean (L_dem x S_dem)
%      mask simply computed from xy and direction.
%  TIPS
%   cmmdl should be defined on the xyz coordinate (x: north, y: east, 
%   z: pointing down)


L_dem = MSLDEMdata.hdr.lines; S_dem = MSLDEMdata.hdr.samples;
dx = MSLDEMdata.hdr.map_info.dx; dy = MSLDEMdata.hdr.map_info.dy;
L_im = camim_size(1); S_im = camim_size(2);
cmmdl.get_image_plane_unit_vectors();

dem_northing   = reshape(MSLDEMdata.hdr.y,L_dem,1);
dem_easting    = reshape(MSLDEMdata.hdr.x,S_dem,1);

dem_imFOV_mask_xyd = false(S_dem,L_dem);
right_dir = false(S_dem,L_dem);

% Placeholder for line by line processing.
deml_g = zeros(S_dem,3);
deml_g(:,2) = dem_easting;

% deml_elev = nan(1,S_dem);
MSLDEMdata.fopen_img(); typeName = 'single';
machine = 'ieee-le';
demlp1_elev = fread(MSLDEMdata.fid_img,S_dem,typeName,0,machine);
demlp1_elev(demlp1_elev==MSLDEMdata.hdr.data_ignore_value) = nan;
xyz_resol_p1_diff = diff(demlp1_elev);
xyz_resol_p1_l = nan(S_dem,1); 
xyz_resol_p1_l(1:S_dem-1) = xyz_resol_p1_diff;
xyz_resol_p1_r = nan(S_dem,1);
xyz_resol_p1_r(2:S_dem) = xyz_resol_p1_diff;

xyz_resol_p1 = nan([S_dem,1]);

H = cmmdl.H';
V = cmmdl.V';
A = cmmdl.A';

Hdash_resol_dxdy = 1.05*(abs(cmmdl.Hdash(1)*dx)+abs(cmmdl.Hdash(2)*dy));
Vdash_resol_dxdy = 1.05*(abs(cmmdl.Vdash(1)*dx)+abs(cmmdl.Vdash(2)*dy));
Hdash3_c = 1.05*cmmdl.Hdash(3);
Vdash3_c = 1.05*cmmdl.Vdash(3);

MSLDEMdata.fopen_img();

for l = 1:L_dem
    %==============================================================
    % Read elevation of the lth line from MSLDEMdata
    %==============================================================
    % demlm1_elev = deml_elev;
    deml_elev = demlp1_elev;
    if l~=L_dem
        demlp1_elev = fread(MSLDEMdata.fid_img,S_dem,typeName,0,machine);
        % demlp1_elev = MSLDEMdata.lazyEnviReadl(l+1,0)';
        demlp1_elev(demlp1_elev==MSLDEMdata.hdr.data_ignore_value) = nan;
    else
        demlp1_elev = nan(S_dem,1);
    end
    
    % converting to xyz
    deml_g(:,1) = dem_northing(l);
    deml_g(:,3) = -deml_elev;
    
    %==============================================================
    % Converting to image coordinate.
    %==============================================================
    % deml_g_vec2d_mC = deml_g' - repmat(cmmdl.C',[1,S_dem]);
    deml_g_vec2d_mC = deml_g - cmmdl.C;
    dst2plane_inprod_vec = deml_g_vec2d_mC*A;
    dem_imx_vec2d = (deml_g_vec2d_mC*H) ./ dst2plane_inprod_vec;
    dem_imy_vec2d = (deml_g_vec2d_mC*V) ./ dst2plane_inprod_vec;
    % imgnd_imxy_vec2d = (HV * deml_g_vec2d_mC) ./ dst2plane_inprod_vec;

    % get resolution
    xyz_resol_m1 = xyz_resol_p1;
    xyz_resol_c1 = abs(xyz_resol_p1_l);
    xyz_resol_c2 = abs(xyz_resol_p1_r);
    
    xyz_resol_p1_diff = diff(demlp1_elev);
    % xyz_resol_p1_l = [xyz_resol_p1_diff nan];
    xyz_resol_p1_l(1:S_dem-1) = xyz_resol_p1_diff;
    xyz_resol_p1_r(2:S_dem) = xyz_resol_p1_diff;
    
    xyz_resol_p11 = deml_elev - demlp1_elev;
    xyz_resol_p12 = abs(xyz_resol_p11 + xyz_resol_p1_l);
    xyz_resol_p13 = abs(xyz_resol_p11 - xyz_resol_p1_r);
    xyz_resol_p11 = abs(xyz_resol_p11);
    
    xyz_resol_p1 = max(max(xyz_resol_p11,xyz_resol_p12),xyz_resol_p13);
    xyz_resol_c = max(xyz_resol_c1,xyz_resol_c2);
    
    deml_elevl_resol = max(max(xyz_resol_p1,xyz_resol_c),xyz_resol_m1);
    
    
    %deml_elevl3 = cat(1,demlm1_elev,deml_elev,demlp1_elev);
    %[deml_elevl3_resol] = get_instantaneous_resolution_for_MSLDEM(deml_elevl3);
    resol_Hdash_vec2d = Hdash_resol_dxdy + Hdash3_c.*deml_elevl_resol;
    resol_Vdash_vec2d = Vdash_resol_dxdy + Vdash3_c.*deml_elevl_resol;
    % resol_Vdash_vec2d = 1.05 + (1.05*HVdash3).*deml_elevl_resol;
    % resol_HVdash_vec2d = repmat(1.05.*HVdash_resol_xy,[1,S_dem]) + (1.05*HVdash3).*deml_elevl_resol;

    % get margin
    hs_d_vec2d = cmmdl.hs.*resol_Hdash_vec2d./abs(dst2plane_inprod_vec);
    vs_d_vec2d = cmmdl.vs.*resol_Vdash_vec2d./abs(dst2plane_inprod_vec);
    % hsvs_d_vec2d = hsvs.*resol_HVdash_vec2d./abs(dst2plane_inprod_vec);

    % test
    dem_imFOV_mask_x_vec = and( (-0.5-hs_d_vec2d) < dem_imx_vec2d, dem_imx_vec2d < (S_im-0.5+hs_d_vec2d) );
    dem_imFOV_mask_y_vec = and( (-0.5-vs_d_vec2d) < dem_imy_vec2d, dem_imy_vec2d < (L_im-0.5+vs_d_vec2d) );
    dem_imFOV_mask_xy_vec = and(dem_imFOV_mask_x_vec,dem_imFOV_mask_y_vec);
%     imgnd_FOV_mask_xy_vec = all(...
%         and( (-0.5-hsvs_d_vec2d) < imgnd_imxy_vec2d,...
%              imgnd_imxy_vec2d < ([S_im;L_im]-0.5+hsvs_d_vec2d) )...
%         ,1);

    right_dir_l = (dst2plane_inprod_vec > 0);

    dem_imFOV_mask_xyd_l = and(right_dir_l,dem_imFOV_mask_xy_vec);
    
    dem_imFOV_mask_xyd(:,l) = dem_imFOV_mask_xyd_l;
    right_dir(:,l) = right_dir_l;
    
end

MSLDEMdata.fclose_img();

% imgnd_FOV_mask_xyd = imgnd_FOV_mask_xyd';
% right_dir = right_dir';

% safeguarding for close distance (some pixels located in the
% opposite direction also needs 

safeguard_mask = false(S_dem,L_dem);
for l = 1:L_dem
    if l==1
        imgnd_FOV_mask_xyd_lm1 = false(S_dem,1);
        dem_imFOV_mask_xyd_l = or(...
        or([false; dem_imFOV_mask_xyd(1:S_dem-1,1)] ...
           , dem_imFOV_mask_xyd(:,1) ...
          ) ...
        , [dem_imFOV_mask_xyd(2:S_dem,1); false]);
    else
        dem_imFOV_mask_xyd_l   = imgnd_FOV_mask_xyd_lm1;
        imgnd_FOV_mask_xyd_lm1 = imgnd_FOV_mask_xyd_lp1;
    end
    if l==L_dem
        imgnd_FOV_mask_xyd_lp1 = false(S_dem,1);
    else
        imgnd_FOV_mask_xyd_lp1 = or(...
            or(dem_imFOV_mask_xyd(:,l+1) ...
               , [false; dem_imFOV_mask_xyd(1:S_dem-1,l+1)]) ...
            ,[dem_imFOV_mask_xyd(2:S_dem,l+1); false]);
    end
    safeguard_mask(:,l) = and(~right_dir(:,l),...
        or(dem_imFOV_mask_xyd_l ...
           , or(imgnd_FOV_mask_xyd_lm1,imgnd_FOV_mask_xyd_lp1) ...
           )...
        );
        
end

% [safeguard_mask] = get_safeguard_imFOV_mask_closedistance(...
%      imgnd_FOV_mask_xyd,right_dir);

dem_imFOV_mask = or(dem_imFOV_mask_xyd,safeguard_mask);

dem_imFOV_mask_xyd = dem_imFOV_mask_xyd';
dem_imFOV_mask = dem_imFOV_mask';


end