function [DDRprj] = proj_crismDDR2mastcam(crism_DDRdata,mastcamdata_obj,varargin)
% [ddr_imFOV_mask] = get_crismDDR_imxyFOV(crism_DDRdata,rover_nav_coord,cmmdl,im_size,varargin)
%   evaluate FOV of an image on an ortho-georeferenced image using a
%   georeferenced DEM image.
%  INPUTS:
%    crism_DDRdata: crism DDRdata
%    mastcamdata_obj: MASTCAMdata or MASTCAMgroup_eye class object
%  OUTPUTS
%   DDRprj: struct having fields:
%    imxy: [L_ddr x S_ddr x 2] 
%      page 1 is the x coordinate in the image frame,
%      page 2 is the y coordinate in the image frame.
%      x and y coordinate values outside of the range between -200 and 
%      200+image_edge are replaced with nans because the values outside of
%      this is far away from the calibration range, therefore accuracy is
%      not guaranteed.
%    imFOV_mask: boolean image, [L_ddr x S_ddr x 1]
%      true if in the FOV, false otherwise.
%      FOV should be true if the linear transformation shows the value
%      between -200 and 200+image_edge. The values outside of this are
%      considered to be inaccurate or invalid. In addition, pixels within
%      the distance 50m are considered to be in FOV, just in case.
%    NorthEastElevation: [L_ddr x S_ddr x 3] 
%      pages 1,2,3 are coordinate in north-east-elevation. 

is_gpu = false;
precision = 'double';
proc_mode = 'BATCH';
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'GPU'
                is_gpu = varargin{i+1};
            case 'PRECISION'
                precision = varargin{i+1};
            case 'PROC_MODE'
                proc_mode = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

if is_gpu
    gpu_varargin = {'gpuArray'};
else
    gpu_varargin = {};
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

if is_gpu
    cmmdl_C = gpuArray(cmmdl_C);
    cmmdl_A = gpuArray(cmmdl_A);
    cmmdl_H = gpuArray(cmmdl_H);
    cmmdl_V = gpuArray(cmmdl_V);
end
cmmdl_HV_mat = cat(1,cmmdl_H,cmmdl_V);

L_ddr = crism_DDRdata.hdr.lines; S_ddr = crism_DDRdata.hdr.samples;
if isempty(crism_DDRdata.ddr), crism_DDRdata.readimg(); end
ddr_geo_latitude_pc = crism_DDRdata.ddr.Latitude.img;
ddr_geo_longitude = crism_DDRdata.ddr.Longitude.img;
ddr_geo_elevation = crism_DDRdata.ddr.Elevation.img;
% planetocentric coordinate to northing easting coordinates
Re = 3396190; % meters ellipsoid radius
ddr_geo_northing = Re .* (pi/180) .* ddr_geo_latitude_pc;
ddr_geo_easting  = Re .* (pi/180) .* ddr_geo_longitude;

if is_gpu
    rov_rot_mat = gpuArray(rov_rot_mat);
    ddr_geo_northing = gpuArray(ddr_geo_northing);
    ddr_geo_easting = gpuArray(ddr_geo_easting);
    ddr_geo_elevation = gpuArray(ddr_geo_elevation);
end

switch lower(precision)
    case 'double'
        ddr_geo_northing = double(ddr_geo_northing);
        ddr_geo_easting = double(ddr_geo_easting);
        ddr_geo_elevation = double(ddr_geo_elevation);
    case 'single'
        ddr_geo_northing = single(ddr_geo_northing);
        ddr_geo_easting = single(ddr_geo_easting);
        ddr_geo_elevation = single(ddr_geo_elevation);
end

ddr_rov0_northing = ddr_geo_northing - rov_northing;
ddr_rov0_easting = ddr_geo_easting - rov_easting;
ddr_rov0_z =  -ddr_geo_elevation - (-rov_elevation); % z is looking down, zenith.

margin_default = 200;

switch upper(proc_mode)
    case 'LINEBYLINE'
        error('Need to check the implementation in this LINEBYLINEMODE');
%         ddr_rov0 = permute(cat(3,ddr_rov0_northing,ddr_rov0_easting,ddr_rov0_z),[3,2,1]);
%         % Perform line by line operation
%         ddr_imxy = nan(L_ddr,S_ddr,2);
%         ddr_imFOV_mask = false(L_ddr,S_ddr,gpu_varargin{:});
%         for l = 1:L_ddr
%             %==========================================================================
%             % projection of ground reference coordinates to ROVER_NAV coordinate
%             %==========================================================================
%             ddrl_rov0 = ddr_rov0(:,:,l);
% 
%             % ---------------------------------------------------------------------
%             % rover coordinate (performing rotation)
%             ddrl_rov = rov_rot_mat \ ddrl_rov0;  % [3 x S_ddr]
% 
%             %==========================================================================
%             % Projection of ROVER_NAV to the image xy coordinate
%             %==========================================================================
%             ddrl_rov_m_cmmdl_C = ddrl_rov-cmmdl_C;
%             right_dir = (cmmdl_A * ddrl_rov_m_cmmdl_C) > -1; % safeguard
% 
%             ddrl_imxy = (cmmdl_HV_mat * ddrl_rov_m_cmmdl_C) ./ (cmmdl_A * ddrl_rov_m_cmmdl_C); % 2 x S_ddr
%             ddrl_imFOV = and(and(all(ddrl_imxy>-200,1),ddrl_imxy(1,:)<S_im+200),ddrl_imxy(2,:)<L_im+200);
%             ddrl_imFOV_mask = and(right_dir,ddrl_imFOV);
% 
%             ddrl_imxy(1,ddrl_imFOV_mask==0) = nan;
%             ddrl_imxy(2,ddrl_imFOV_mask==0) = nan;
% 
%             ddr_imxy(l,:,:) = permute(ddrl_imxy, [3,2,1]);
% 
%             ddr_imFOV_mask(l,:) = ddrl_imFOV_mask;
%         end
    case 'BATCH'
        ddr_rov0 = permute(cat(3,ddr_rov0_northing,ddr_rov0_easting,ddr_rov0_z),[3,1,2]);
        ddr_rov0_2d = reshape(ddr_rov0,[3,L_ddr*S_ddr]);
        % rover coordinate (performing rotation)
        ddr_rov_2d = rov_rot_mat \ ddr_rov0_2d;
        ddr_rov_m_cmmdl_C = ddr_rov_2d-cmmdl_C;
        

        ddr_imxy = (cmmdl_HV_mat * ddr_rov_m_cmmdl_C) ./ (cmmdl_A * ddr_rov_m_cmmdl_C); % 2 x S_ddr
        marg = margin_default;
        % sg_mag = 1;
        ddr_imFOV = and(and(all(ddr_imxy>-marg,1),ddr_imxy(1,:)<S_im+marg),...
                ddr_imxy(2,:)<L_im+marg);
        
        right_dir = (cmmdl_A * ddr_rov_m_cmmdl_C) > 0; %(-1)*sg_mag; % safeguard
        
        ddr_imFOV_mask = and(right_dir,ddr_imFOV);
        
        ddr_imxy(1,ddr_imFOV_mask==0) = nan;
        ddr_imxy(2,ddr_imFOV_mask==0) = nan;
        
        ddr_imFOV_tooclose = sqrt(sum(ddr_rov_m_cmmdl_C.^2,1)) < 100;
        ddr_imFOV_mask = or(ddr_imFOV_mask,ddr_imFOV_tooclose);       
        ddr_imFOV_mask = reshape(ddr_imFOV_mask,[L_ddr,S_ddr]);
         
    otherwise
        error('undefined processing mode: %s',proc_mode);
        
end

% safeguard
for l=1:L_ddr
    ddr_imFOV_mask(l,:) = conv(ddr_imFOV_mask(l,:),[1 1 1],'same')>0;
end

for s=1:S_ddr
    ddr_imFOV_mask(:,s) = conv(ddr_imFOV_mask(:,s),[1;1;1],'same')>0;
end

% ddr_imxy(1,ddr_imFOV_mask==0) = nan;
% ddr_imxy(2,ddr_imFOV_mask==0) = nan;
ddr_imxy = permute(reshape(ddr_imxy,[2,L_ddr,S_ddr]),[2,3,1]);
ddr_geo = cat(3,ddr_geo_northing,ddr_geo_easting,ddr_geo_elevation);

DDRprj = [];
DDRprj.imxy = ddr_imxy;
DDRprj.NorthEastElevation = ddr_geo;
DDRprj.imFOV_mask = ddr_imFOV_mask;


end