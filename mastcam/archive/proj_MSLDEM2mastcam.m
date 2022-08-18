function [MSLDEMprj] = proj_MSLDEM2mastcam(MSLDEMdata,mastcamdata_obj,varargin)
% [MSLDEMprj] = get_dem_imFOV(basename_dem,dpath_dem,rover_nav_coord,cmmdl,im_size,varargin)
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
is_gpu = 0; precision = 'double';
PROC_MODE = 'naive';
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'PROC_MODE'
                PROC_MODE = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

%% GET CAMERA, ROVER_NAV, and Image size information
%-------------------------------------------------------------------------%
% Get camera CAHV model parameters
%-------------------------------------------------------------------------%

cmmdl = mastcamdata_obj.CAM_MDL;
cmmdl_C = reshape(cmmdl.C,3,1);
cmmdl_A = reshape(cmmdl.A,1,3);
cmmdl_H = reshape(cmmdl.H,1,3);
cmmdl_V = reshape(cmmdl.V,1,3);

if is_gpu
    cmmdl_C = gpuArray(cmmdl_C);
    cmmdl_A = gpuArray(cmmdl_A);
    cmmdl_H = gpuArray(cmmdl_H);
    cmmdl_V = gpuArray(cmmdl_V);
end

cmmdl_HV_mat = cat(1,cmmdl_H,cmmdl_V);

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
% Get Rover Nav
%-------------------------------------------------------------------------%
rover_nav_coord = mastcamdata_obj.ROVER_NAV;
rov_northing  = rover_nav_coord.NORTHING;
rov_easting   = rover_nav_coord.EASTING;
rov_elevation = rover_nav_coord.ELEVATION;
rov_rot_mat   = rover_nav_coord.rot_mat;

%% Get MSL DEM information
%-------------------------------------------------------------------------%
% Get MSL DEM information
%-------------------------------------------------------------------------%
L_dem = MSLDEMdata.hdr.lines; S_dem = MSLDEMdata.hdr.samples;
dem_northing   = reshape(MSLDEMdata.hdr.y,L_dem,1);
dem_easting    = reshape(MSLDEMdata.hdr.x,1,S_dem);

%-------------------------------------------------------------------------%
% Northing easting to rover coordinate with no rotation. Pre-computation
% outside of the for-loops
%-------------------------------------------------------------------------%
dem_rov0_x = dem_northing - rov_northing;
dem_rov0_y = dem_easting - rov_easting;

%% First compute dem_imFOV_mask
switch lower(PROC_MODE)
    case 'naive'
        marg = 300;
        %-----------------------------------------------------------------%
        % Placeholders for processing
        %-----------------------------------------------------------------%
        dem_imFOV_mask_cmmdl = false(L_dem,S_dem);
        % dem_imFOV_mask_cmmdl: 
        %   Solely computed from camera model. If the coordinate values is 
        %   too away from the field of view, then it is considered as 
        %   invalid.

        dem_imFOV_mask = false(L_dem,S_dem);
        % dem_imFOV_mask: 
        %   Some safe guard is applied. Close pixels are forcefully set in 
        %   the FOV as such pixels are easily get too far away in the 
        %   coordinate, but these pixels may be necessary for calculating 
        %   the intersection.

        % Placeholder for line by line processing.
        deml_rov0 = zeros(3,S_dem,precision);
        deml_rov0(2,:) = dem_rov0_y;
        
        %-----------------------------------------------------------------%
        % Perform line by line operation
        for l = 1:L_dem
            %==============================================================
            % Read elevation of the lth line from MSLDEMdata
            %==============================================================
            deml_elevation = MSLDEMdata.lazyEnviReadl(l,0);
            deml_elevation(deml_elevation==MSLDEMdata.hdr.data_ignore_value) = nan;

            %==============================================================
            % projection of ground reference coordinates to ROVER_NAV coordinate
            %==============================================================
            % northing easting to rover coordinate with no rotation
            deml_rov0(1,:) = dem_rov0_x(l);
            deml_rov0(3,:) = -deml_elevation - (-rov_elevation);

            % rover coordinate (performing rotation)
            deml_rov = rov_rot_mat \ deml_rov0;  % [3 x S_geo]

            %==============================================================
            % Projection of ROVER_NAV to the image xy coordinate
            %==============================================================
            deml_rov_m_cmmdl_C = deml_rov-cmmdl_C;
            deml_im = (cmmdl_HV_mat * deml_rov_m_cmmdl_C) ./ (cmmdl_A * deml_rov_m_cmmdl_C);
            % 2 x S_geo
            %==============================================================
            % Evaluate Field of View (FOV)
            %==============================================================
            % Some safeguarding with margin. Valid values are valid image
            % coordinate values plus or minus marg (default 200).
            deml_imFOV = and(and(all(deml_im>-marg,1),deml_im(1,:)<S_im+marg),...
                deml_im(2,:)<L_im+marg);

            % make sure diretion is right in the direction of axis vector
            right_dir = (cmmdl_A * deml_rov_m_cmmdl_C) > -1;
            deml_imFOV_mask = and(right_dir,deml_imFOV);

            % If the pixel close to the camera center might need to be complementd.
            deml_imFOV_tooclose = sqrt(sum(deml_rov_m_cmmdl_C.^2,1)) < 5;
            deml_imFOV_mask_comb = or(deml_imFOV_mask,deml_imFOV_tooclose);

            dem_imFOV_mask_cmmdl(l,:) = deml_imFOV_mask;
            dem_imFOV_mask(l,:) = deml_imFOV_mask_comb;
        end
        
        % Final safeguarding. Make sure the FOV is completely covered by the mask
        % by expanding the mask one pixel more.
        for l=1:L_dem
            dem_imFOV_mask(l,:) = conv(dem_imFOV_mask(l,:),[1 1 1],'same')>0;
        end
        for s=1:S_dem
            dem_imFOV_mask(:,s) = conv(dem_imFOV_mask(:,s),[1;1;1],'same')>0;
        end
        
    case 'tan3'
        % convert camera model in the ground reference coordinate
        
        %-----------------------------------------------------------------%
        % Get some camera parameters used for later processing.
        %-----------------------------------------------------------------%
        % get FOV in the horizontal and vertical directions
        % probably the the most upper left pixel is defined as (0 0)
        % [fovh,fovv] = cmmdl.get_FOV([0 S_im]-0.5,[0 L_im]-0.5);
        % tan_fovh = (tan(fovh(1)) + tan(fovh(2)));
        % tan_fovv = (tan(fovv(1)) + tan(fovv(2)));
        % get unit vectors in the x and y directions in the ROVER_NAV
        % coordinate.
        cmmdl.get_image_plane_unit_vectors();
        % convert them in the north-east-elevation coordinate.
        Hdash_rov0 = cmmdl.Hdash * rov_rot_mat';
        Vdash_rov0 = cmmdl.Vdash * rov_rot_mat';
        H_rov0 = cmmdl.H * rov_rot_mat';
        V_rov0 = cmmdl.V * rov_rot_mat';
        A_rov0 = cmmdl.A * rov_rot_mat';
        
        
        %-----------------------------------------------------------------%
        % Placeholders for processing
        %-----------------------------------------------------------------%
        dem_imFOV_mask_cmmdl = false(L_dem,S_dem);
        % dem_imFOV_mask_cmmdl: 
        %   Solely computed from camera model. If the coordinate values is 
        %   too away from the field of view, then it is considered as 
        %   invalid.
        
        right_dir = false(L_dem,S_dem);

        % Placeholder for line by line processing.
        deml_rov0 = zeros(3,S_dem,precision);
        deml_rov0(2,:) = dem_rov0_y;
        
        %-----------------------------------------------------------------%
        % Perform line by line operation
        deml_elevation = nan(1,S_dem);
        demlp1_elevation = MSLDEMdata.lazyEnviReadl(1,0);
        demlp1_elevation(demlp1_elevation==MSLDEMdata.hdr.data_ignore_value) = nan;
        for l = 1:L_dem
            %==============================================================
            % Read elevation of the lth line from MSLDEMdata
            %==============================================================
            demlm1_elevation = deml_elevation;
            deml_elevation = demlp1_elevation;
            demlp1_elevation = MSLDEMdata.lazyEnviReadl(l+1,0);
            demlp1_elevation(demlp1_elevation==MSLDEMdata.hdr.data_ignore_value) = nan;
            
            %==============================================================
            % projection of ground reference coordinates to ROVER_NAV coordinate
            %==============================================================
            % northing easting to rover coordinate with no rotation
            deml_rov0(1,:) = dem_rov0_x(l);
            deml_rov0(3,:) = -deml_elevation - (-rov_elevation);

            % rover coordinate (performing rotation)
            deml_rov = rov_rot_mat \ deml_rov0;  % [3 x S_dem]
            
            %==============================================================
            % Projection of ROVER_NAV to the image xy coordinate
            %==============================================================
            deml_rov_m_cmmdl_C = deml_rov-cmmdl_C;
            deml_imxy = (cmmdl_HV_mat * deml_rov_m_cmmdl_C) ./ (cmmdl_A * deml_rov_m_cmmdl_C);
            
            dst2plane_inprod = (cmmdl_A * deml_rov_m_cmmdl_C);
            dst2plane = abs(dst2plane_inprod);
            wdth = dst2plane .* tan_fovh;
            hgt  = dst2plane .* tan_fovv;
            deml_elevation_resol = get_instantaneous_resolution(...
                cat(2,demlm1_elevation,deml_elevation,demlp1_elevation) );
            
            mrgn_x = S_im./wdth.*(...
                abs(Hdash_rov0(1)) * MSLDEMdata.hdr.map_info.dx ...
                + abs(Hdash_rov0(2)) * MSLDEMdata.hdr.map_info.dy ...
                + abs(Hdash_rov0(3)) * deml_elevation_resol(2,:));
            mrgn_y = L_im./hgt.*(...
                abs(Vdash_rov0(1)) * MSLDEMdata.hdr.map_info.dx ...
                + abs(Vdash_rov0(2)) * MSLDEMdata.hdr.map_info.dy ...
                + abs(Vdash_rov0(3)) * deml_elevation_resol(2,:));

            %==============================================================
            % Evaluate Field of View (FOV)
            %==============================================================
            deml_imFOV_mask_x = and(-mrgn_x<deml_imxy(1,:),deml_imxy(1,:)<S_im+mrgn_x);
            deml_imFOV_mask_y = and(-mrgn_y<deml_imxy(2,:),deml_imxy(2,:)<S_im+mrgn_y);
            deml_imFOV_mask_xy = and(deml_imFOV_mask_x,deml_imFOV_mask_y);
            
            right_dir_l = (dst2plane_inprod > 0);
            deml_imFOV_mask_xyd = and(right_dir_l,deml_imFOV_mask_xy);

            dem_imFOV_mask_cmmdl(l,:) = deml_imFOV_mask_xyd;
            right_dir(l,:) = right_dir_l;
            
        end
        
        % safeguarding for close distance (some pixels located in the
        % opposite direction also needs 
        dem_imFOV_mask_cmmdl_pd = padarray(dem_imFOV_mask_cmmdl,[1 1],'replicate','both');

        dem_imFOV_mask_cmmdl_nbr = false(L_dem,S_dem);
        for i=1:3
            for j=1:3
                dem_imFOV_mask_cmmdl_nbr = or(dem_imFOV_mask_cmmdl_nbr,dem_imFOV_mask_cmmdl_pd(i:end-3+i,j:end-3+j));
            end
        end
        safeguard_mask = and(~right_dir,dem_imFOV_mask_cmmdl_nbr);
        dem_imFOV_mask = or(dem_imFOV_mask_cmmdl,safeguard_mask);
        
        clear safeguard_mask right_dir dem_imFOV_mask_cmmdl_pd dem_imFOV_mask_cmmdl_nbr
end



%% Next computing image coordinate for each pixel of DEM image.
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
dem_imxy = nan(len_vl,len_vs,2);

dem_imFOV_mask_cmmdl_vls = dem_imFOV_mask_cmmdl(line_range,sample_range);
% pre-asign easting for each line of the data to be analyzed.
deml_rov0_vs = zeros(3,len_vs,precision);
deml_rov0_vs(2,:) = dem_rov0_y(sample_range);
for li = 1:len_vl
    l = line_range(li);
    
    %======================================================================
    % Read elevation of the lth line from MSLDEMdata
    %======================================================================
    deml_elevation = MSLDEMdata.lazyEnviReadl(l,0);
    deml_elevation(deml_elevation==MSLDEMdata.hdr.data_ignore_value) = nan;
    
    %======================================================================
    % projection of ground reference coordinates to ROVER_NAV coordinate
    %======================================================================
    % northing easting to rover coordinate with no rotation
    deml_rov0_vs(1,:) = dem_rov0_x(l);
    deml_rov0_vs(3,:) = -deml_elevation(sample_range) - (-rov_elevation);
    
    % rover coordinate (performing rotation)
    deml_rov_vs = rov_rot_mat \ deml_rov0_vs;  % [3 x S_geo]
    
    %======================================================================
    % Projection of ROVER_NAV to the image xy coordinate
    %======================================================================
    deml_rov_vs_m_cmmdl_C = deml_rov_vs-cmmdl_C;
    deml_im = (cmmdl_HV_mat * deml_rov_vs_m_cmmdl_C) ./ ...
        (cmmdl_A * deml_rov_vs_m_cmmdl_C); % 2 x S_geo
    
    % Some safeguarding with margin. Valid values are valid image
    % coordinate values plus or minus marg (default 200). So if the values
    % are out of this range, then NaNs are substituted.
    deml_im(1,dem_imFOV_mask_cmmdl_vls(li,:)==0) = nan;
    deml_im(2,dem_imFOV_mask_cmmdl_vls(li,:)==0) = nan;
    
    dem_imxy(li,:,:) = permute(deml_im, [3,2,1]);
    
end

hdr_dem_imxy = MSLDEMdata.hdr;
hdr_dem_imxy.interleave = 'bil';
hdr_dem_imxy.bands = 2;
hdr_dem_imxy.lines = len_vl;
hdr_dem_imxy.samples = len_vs;
hdr_dem_imxy.band_names = {'imx','imy'};
hdr_dem_imxy.line_offset = line_offset;
hdr_dem_imxy.sample_offset = sample_offset;

%% Summary
MSLDEMdata.fclose_img();

MSLDEMprj = [];
MSLDEMprj.imFOV_mask = dem_imFOV_mask;
MSLDEMprj.imxy = dem_imxy;
MSLDEMprj.hdr_imxy = hdr_dem_imxy;

end