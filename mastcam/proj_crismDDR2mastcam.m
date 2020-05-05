function [DDRprj] = proj_crismDDR2mastcam(crism_DDRdata,mastcamdata_obj,varargin)
% [DDRprj] = proj_crismDDR2mastcam(crism_DDRdata,mastcamdata_obj,varargin)
%   evaluate FOV of an MASTCAM image on the CRISM image and calculate
%   corresponding (x,y) coordinates of pixels of CRISM image in the MASTCAM
%   image coordinate system.
%  INPUTS:
%    crism_DDRdata: crism DDRdata
%    mastcamdata_obj: MASTCAMdata or MASTCAMgroup_eye class object
%  OUTPUTS
%   DDRprj: struct having fields:
%    imxy: [L_ddr x S_ddr x 2] 
%      page 1 is the x coordinate in the image frame,
%      page 2 is the y coordinate in the image frame. 
%    imFOV_mask: boolean image, [L_ddr x S_ddr x 1]
%      true if in the FOV, false otherwise.
%      FOV should be true if the linear transformation shows the value
%      between -200 and 200+image_edge. The values outside of this are
%      considered to be inaccurate or invalid. In addition, pixels within
%      the distance 50m are considered to be in FOV, just in case.
%    NorthEastElevation: [L_ddr x S_ddr x 3] 
%      pages 1,2,3 are coordinate in north-east-elevation. 
%  OPTIONAL PARAMETERS
%    'PROC_MODE': {'naive','tan','tan2'}
%      - 'naive -
%      x and y coordinate values outside of the range between -[MARGIN] and 
%      [MARGIN]+image_edge are replaced with nans because the values 
%      outside of this is far away from the calibration range, therefore 
%      accuracy is not guaranteed. [MARGIN] is a constant value and set to 
%      200 by default and can be changed by setting 'MARGIN' parameter. 
%      In addition, pixels close to the camera position is also added.
%      
%      - 'tan' -
%      x and y coordinate values outside of the range between -[MARGIN] and 
%      [MARGIN]+image_edge are replaced with nans. 
%
%      - 'tan2' -
%      x and y coordinate values outside of the range between -[MARGIN] and 
%      [MARGIN]+image_edge are replaced with nans.
%      

mrgn_default = 200;
PROC_MODE = 'naive';
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'MARGIN'
                mrgn_default = varargin{i+1};
            case 'PROC_MODE'
                PROC_MODE = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

cmmdl = mastcamdata_obj.CAM_MDL;
cmmdl_C = reshape(cmmdl.C,3,1);
cmmdl_A = reshape(cmmdl.A,1,3);
cmmdl_H = reshape(cmmdl.H,1,3);
cmmdl_V = reshape(cmmdl.V,1,3);

% Get rover coordinate
rover_nav_coord = mastcamdata_obj.ROVER_NAV;
rov_northing  = rover_nav_coord.NORTHING;
rov_easting   = rover_nav_coord.EASTING;
rov_elevation = rover_nav_coord.ELEVATION;
rov_rot_mat   = rover_nav_coord.rot_mat;

switch class(mastcamdata_obj)
    case 'MASTCAMdata'
        L_im = mastcamdata_obj.hdr.lines; S_im = mastcamdata_obj.samples;
    case 'MASTCAMgroup_eye'
        L_im = mastcamdata_obj.L_im; S_im = mastcamdata_obj.S_im;
    otherwise
        error('second input needs to be either MASTCAMdata or MASTCAMgroup_eye class.');
end

cmmdl_HV_mat = cat(1,cmmdl_H,cmmdl_V);

L_ddr = crism_DDRdata.hdr.lines; S_ddr = crism_DDRdata.hdr.samples;
[ddr_nee] = get_CRISM_NorthEastElevation(crism_DDRdata,'RE',3396190);

%% MAIN processing
%==========================================================================
% Projection of ground reference coordinates to ROVER_NAV coordinate
%==========================================================================
ddr_rov0 = ddr_nee;
ddr_rov0(:,:,1) = ddr_rov0(:,:,1) - rov_northing;
ddr_rov0(:,:,2) = ddr_rov0(:,:,2) - rov_easting;
ddr_rov0(:,:,3) = -ddr_rov0(:,:,3) - (-rov_elevation); % z is looking down, zenith.
% ddr_rov0_northing = ddr_nee(:,:,1) - rov_northing;
% ddr_rov0_easting = ddr_nee(:,:,2) - rov_easting;
% ddr_rov0_z =  -ddr_nee(:,:,3) - (-rov_elevation); % z is looking down, zenith.
% ddr_rov0 = permute(cat(3,ddr_rov0_northing,ddr_rov0_easting,ddr_rov0_z),[3,1,2]);
ddr_rov0_vec2d = reshape(ddr_rov0,[L_ddr*S_ddr,3])'; % 3 x ( L_ddr*S_ddr )
% rover coordinate (performing rotation)
ddr_rov_vec2d = rov_rot_mat \ ddr_rov0_vec2d;

%==========================================================================
% Projection of ROVER_NAV to the image xy coordinate
%==========================================================================
ddr_rov_m_cmmdl_C_vec2d = ddr_rov_vec2d-cmmdl_C; % 3 x ( L_ddr*S_ddr )
dst2plane_inprod_vec = cmmdl_A*ddr_rov_m_cmmdl_C_vec2d;
ddr_imxy_vec2d = (cmmdl_HV_mat * ddr_rov_m_cmmdl_C_vec2d)./dst2plane_inprod_vec;
% ddr_imxy_vec2d: 2 x ( L_ddr*S_ddr )


switch lower(PROC_MODE)
    case 'naive'
        mrgn = mrgn_default;
        % sg_mag = 1;
        ddr_imFOV_mask_xy_vec = and(and(all(ddr_imxy_vec2d>-mrgn,1),ddr_imxy_vec2d(1,:)<S_im+mrgn),...
                ddr_imxy_vec2d(2,:)<L_im+mrgn);

        right_dir_vec = (cmmdl_A * ddr_rov_m_cmmdl_C_vec2d) > 0; %(-1)*sg_mag; % safeguard

        ddr_imFOV_mask_xyd = and(right_dir_vec,ddr_imFOV_mask_xy_vec);

        ddr_imxy_vec2d(1,ddr_imFOV_mask_xyd==0) = nan;
        ddr_imxy_vec2d(2,ddr_imFOV_mask_xyd==0) = nan;

        ddr_imFOV_tooclose = sqrt(sum(ddr_rov_m_cmmdl_C_vec2d.^2,1)) < 100;
        ddr_imFOV_mask = or(ddr_imFOV_mask_xyd,ddr_imFOV_tooclose);       
        ddr_imFOV_mask = reshape(ddr_imFOV_mask,[L_ddr,S_ddr]);
        
        % safeguard
        for l=1:L_ddr
            ddr_imFOV_mask(l,:) = conv(ddr_imFOV_mask(l,:),[1 1 1],'same')>0;
        end

        for s=1:S_ddr
            ddr_imFOV_mask(:,s) = conv(ddr_imFOV_mask(:,s),[1;1;1],'same')>0;
        end
        
    case 'tan'
        %==================================================================
        % 
        %==================================================================
        % get FOV in the horizontal and vertical directions
        % probably the the most upper left pixel is defined as (0 0)
        % [fovh,fovv] = cmmdl.get_FOV([0 S_im]-0.5,[0 L_im]-0.5);
        % get Distance to the image plane where each pixel lies.
        % dst2plane_inprod_vec = cmmdl_A*ddr_rov_m_cmmdl_C_vec2d;
        dst2plane_vec = abs(dst2plane_inprod_vec);
        % Get the width and hight for each image plane of distance.
        % wdth = dst2plane_vec .* (tan(fovh(1)) + tan(fovh(2)));
        % hgt = dst2plane_vec .* (tan(fovv(1)) + tan(fovv(2)));
        
        cmmdl.get_image_plane_unit_vectors();
        wdth = dst2plane_vec * (S_im/cmmdl.hs);
        hgt = dst2plane_vec * (L_im/cmmdl.vs);

        % Get approximated resolution for each pixel in the (x,y,z)
        % directions.
        [ddr_rov_resol] = get_instantaneous_resolution(...
            reshape(ddr_rov_vec2d',[L_ddr,S_ddr,3]));
        ddr_rov_resol_vec2d = reshape(ddr_rov_resol,[L_ddr*S_ddr,3])';
        % get unit vectors in the x and y directions.
        [Hdash,Vdash] = cmmdl.get_image_plane_unit_vectors();
        
        % Get margins 
        mrgn_x_vec = S_im./wdth.*(abs(Hdash)*ddr_rov_resol_vec2d);
        mrgn_y_vec = L_im./hgt.*(abs(Vdash)*ddr_rov_resol_vec2d);
        
        ddr_imFOV_mask_x_vec = and(...
            -mrgn_x_vec<ddr_imxy_vec2d(1,:),ddr_imxy_vec2d(1,:)<S_im+mrgn_x_vec);
        ddr_imFOV_mask_y_vec = and(...
            -mrgn_y_vec<ddr_imxy_vec2d(2,:),ddr_imxy_vec2d(2,:)<S_im+mrgn_y_vec);
        ddr_imFOV_mask_xy_vec = and(ddr_imFOV_mask_x_vec,ddr_imFOV_mask_y_vec);

        right_dir_vec = (dst2plane_inprod_vec > 0);

        ddr_imFOV_mask_xyd = and(right_dir_vec,ddr_imFOV_mask_xy_vec);
        
        ddr_imxy_vec2d(1,ddr_imFOV_mask_xyd==0) = nan;
        ddr_imxy_vec2d(2,ddr_imFOV_mask_xyd==0) = nan;
        
        % safeguarding for close distance (some pixels located in the
        % opposite direction also needs 
        ddr_imFOV_mask_xyd = reshape(ddr_imFOV_mask_xyd,[L_ddr,S_ddr]);
        ddr_imFOV_mask_xyd_pd = padarray(ddr_imFOV_mask_xyd,[1 1],'replicate','both');
        
        ddr_imFOV_mask_xyd_nbr = false(L_ddr,S_ddr);
        for i=1:3
            for j=1:3
                ddr_imFOV_mask_xyd_nbr = or(ddr_imFOV_mask_xyd_nbr,ddr_imFOV_mask_xyd_pd(i:end-3+i,j:end-3+j));
            end
        end
        safeguard_mask = and(reshape(right_dir_vec==0,L_ddr,S_ddr),ddr_imFOV_mask_xyd_nbr);
        ddr_imFOV_mask = or(ddr_imFOV_mask_xyd,safeguard_mask);
        
    case 'tan2'
        %==================================================================
        % 
        %==================================================================
        dst2plane_vec = abs(dst2plane_inprod_vec);
        right_dir_vec_sgn = sign(dst2plane_inprod_vec);
        % get unit vectors in the x and y directions.
        cmmdl.get_image_plane_unit_vectors();
        % get FOV in the horizontal and vertical directions
        % probably the the most upper left pixel is defined as (0 0)
        
        % [pmc_h,pmc_v] = cmmdl.get_pmc_FOV([0 S_im]-0.5,[0 L_im]-0.5);

        % x_edge = (cmmdl.Hdash*pmc_h./(cmmdl.A * pmc_h))' .* dst2plane_vec;
        % y_edge = (cmmdl.Vdash*pmc_v./(cmmdl.A * pmc_v))' .* dst2plane_vec;
        
        x_edge = ([0;S_im]-0.5-cmmdl.hc)./cmmdl.hs .* dst2plane_vec;
        y_edge = ([0;L_im]-0.5-cmmdl.vc)./cmmdl.vs .* dst2plane_vec;
        
        % Get approximated resolution for each pixel in the (x,y,z)
        % directions.
        [ddr_rov_resol] = get_instantaneous_resolution(...
            reshape(ddr_rov_vec2d',[L_ddr,S_ddr,3]));
        ddr_rov_resol_vec2d = reshape(ddr_rov_resol,[L_ddr*S_ddr,3])';
        
        resol_Hdash = (abs(cmmdl.Hdash)*ddr_rov_resol_vec2d);
        resol_Vdash = (abs(cmmdl.Vdash)*ddr_rov_resol_vec2d);
        
        p_Hdash = (cmmdl.Hdash * ddr_rov_m_cmmdl_C_vec2d) .* right_dir_vec_sgn;
        p_Vdash = (cmmdl.Vdash * ddr_rov_m_cmmdl_C_vec2d) .* right_dir_vec_sgn;
        
        ddr_imFOV_mask_x_vec = and((x_edge(1,:)-resol_Hdash)<p_Hdash,p_Hdash<(x_edge(2,:) + resol_Hdash));
        ddr_imFOV_mask_y_vec = and((y_edge(1,:)-resol_Vdash)<p_Vdash,p_Vdash<(y_edge(2,:) + resol_Vdash));
        ddr_imFOV_mask_xy_vec = and(ddr_imFOV_mask_x_vec,ddr_imFOV_mask_y_vec);

        right_dir_vec = (dst2plane_inprod_vec > 0);

        ddr_imFOV_mask_xyd = and(right_dir_vec,ddr_imFOV_mask_xy_vec);
        
        ddr_imxy_vec2d(1,ddr_imFOV_mask_xyd==0) = nan;
        ddr_imxy_vec2d(2,ddr_imFOV_mask_xyd==0) = nan;
        
        % safeguarding for close distance (some pixels located in the
        % opposite direction also needs 
        ddr_imFOV_mask_xyd = reshape(ddr_imFOV_mask_xyd,[L_ddr,S_ddr]);
        ddr_imFOV_mask_xyd_pd = padarray(ddr_imFOV_mask_xyd,[1 1],'replicate','both');
        
        ddr_imFOV_mask_xyd_nbr = false(L_ddr,S_ddr);
        for i=1:3
            for j=1:3
                ddr_imFOV_mask_xyd_nbr = or(ddr_imFOV_mask_xyd_nbr,ddr_imFOV_mask_xyd_pd(i:end-3+i,j:end-3+j));
            end
        end
        safeguard_mask = and(reshape(right_dir_vec==0,L_ddr,S_ddr),ddr_imFOV_mask_xyd_nbr);
        ddr_imFOV_mask = or(ddr_imFOV_mask_xyd,safeguard_mask);
        
    case 'tan3'
        cmmdl.get_image_plane_unit_vectors();
        dst2plane_vec = abs(dst2plane_inprod_vec);
        
        % get resolution
        [ddr_rov_resol] = get_instantaneous_resolution(...
            reshape(ddr_rov_vec2d',[L_ddr,S_ddr,3]));
        ddr_rov_resol_vec2d = reshape(ddr_rov_resol,[L_ddr*S_ddr,3])';
        resol_Hdash = (abs(cmmdl.Hdash)*ddr_rov_resol_vec2d);
        resol_Vdash = (abs(cmmdl.Vdash)*ddr_rov_resol_vec2d);
        
        % get margin
        hs_d = cmmdl.hs*resol_Hdash./dst2plane_vec;
        vs_d = cmmdl.vs*resol_Vdash./dst2plane_vec;
        
        % test
        ddr_imFOV_mask_x_vec = and((-0.5-hs_d)<ddr_imxy_vec2d(1,:),ddr_imxy_vec2d(1,:)<(S_im-0.5+hs_d));
        ddr_imFOV_mask_y_vec = and((-0.5-vs_d)<ddr_imxy_vec2d(2,:),ddr_imxy_vec2d(2,:)<(L_im-0.5+vs_d));
        ddr_imFOV_mask_xy_vec = and(ddr_imFOV_mask_x_vec,ddr_imFOV_mask_y_vec);
        
        right_dir_vec = (dst2plane_inprod_vec > 0);

        ddr_imFOV_mask_xyd = and(right_dir_vec,ddr_imFOV_mask_xy_vec);
        
        ddr_imxy_vec2d(1,ddr_imFOV_mask_xyd==0) = nan;
        ddr_imxy_vec2d(2,ddr_imFOV_mask_xyd==0) = nan;
        
        % safeguarding for close distance (some pixels located in the
        % opposite direction also needs 
        ddr_imFOV_mask_xyd = reshape(ddr_imFOV_mask_xyd,[L_ddr,S_ddr]);
        ddr_imFOV_mask_xyd_pd = padarray(ddr_imFOV_mask_xyd,[1 1],'replicate','both');
        
        ddr_imFOV_mask_xyd_nbr = false(L_ddr,S_ddr);
        for i=1:3
            for j=1:3
                ddr_imFOV_mask_xyd_nbr = or(ddr_imFOV_mask_xyd_nbr,ddr_imFOV_mask_xyd_pd(i:end-3+i,j:end-3+j));
            end
        end
        safeguard_mask = and(reshape(right_dir_vec==0,L_ddr,S_ddr),ddr_imFOV_mask_xyd_nbr);
        ddr_imFOV_mask = or(ddr_imFOV_mask_xyd,safeguard_mask);
        
        
end
         

%%

DDRprj = [];
DDRprj.imxy = permute(reshape(ddr_imxy_vec2d,[2,L_ddr,S_ddr]),[2,3,1]);
DDRprj.NorthEastElevation = ddr_nee;
DDRprj.imFOV_mask = ddr_imFOV_mask;


end