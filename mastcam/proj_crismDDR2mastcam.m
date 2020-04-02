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

margin_default = 200;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'MARGIN'
                margin_default = varargin{i+1};
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

ddr_rov0_northing = ddr_nee(:,:,1) - rov_northing;
ddr_rov0_easting = ddr_nee(:,:,2) - rov_easting;
ddr_rov0_z =  -ddr_nee(:,:,3) - (-rov_elevation); % z is looking down, zenith.

%% MAIN processing
    
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
         

%%
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
% ddr_geo = cat(3,ddr_geo_northing,ddr_geo_easting,ddr_geo_elevation);

DDRprj = [];
DDRprj.imxy = ddr_imxy;
DDRprj.NorthEastElevation = ddr_nee;
DDRprj.imFOV_mask = ddr_imFOV_mask;


end