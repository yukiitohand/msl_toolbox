function [DDRprj] = proj_crismDDR2mastcam_v2(crism_DDRdata,mastcamdata_obj,varargin)
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
%    'PROC_MODE': {'tan3_geo','tan3_rov'}
%      

PROC_MODE = 'tan3_geo';
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

cmmdl = mastcamdata_obj.CAM_MDL;
rover_nav_coord = mastcamdata_obj.ROVER_NAV;

switch class(mastcamdata_obj)
    case 'MASTCAMdata'
        L_im = mastcamdata_obj.hdr.lines; S_im = mastcamdata_obj.samples;
    case 'MASTCAMgroup_eye'
        L_im = mastcamdata_obj.L_im; S_im = mastcamdata_obj.S_im;
    otherwise
        error('second input needs to be either MASTCAMdata or MASTCAMgroup_eye class.');
end

[ddr_nee] = get_CRISM_NorthEastElevation(crism_DDRdata,'RE',3396190);

%% MAIN processing
switch lower(PROC_MODE)
    case 'tan3_geo'
        ddr_xyz = ddr_nee;
        ddr_xyz(:,:,3) = (-1)*ddr_xyz(:,:,3);
        cmmdl_geo = transform_CAHVOR_MODEL_wROVER_NAV(cmmdl,rover_nav_coord);
        [ddr_imFOV_mask,ddr_imxy] = get_imFOV_mask_imgnd_xyz(cmmdl_geo,[L_im,S_im],ddr_xyz);
        
    case 'tan3_rov'
        L_ddr = crism_DDRdata.hdr.lines; S_ddr = crism_DDRdata.hdr.samples;
        ddr_rov0 = ddr_nee;
        ddr_rov0(:,:,1) = ddr_rov0(:,:,1) - rover_nav_coord.NORTHING;
        ddr_rov0(:,:,2) = ddr_rov0(:,:,2) - rover_nav_coord.EASTING;
        ddr_rov0(:,:,3) = -ddr_rov0(:,:,3) - (-rover_nav_coord.ELEVATION); % z is looking down, zenith.
        % rover coordinate (performing rotation)
        ddr_rov_vec2d = rover_nav_coord.rot_mat_inv * reshape(ddr_rov0,[L_ddr*S_ddr,3])';
        ddr_rov = reshape(ddr_rov_vec2d',[L_ddr,S_ddr,3]);
        
        [ddr_imFOV_mask,ddr_imxy] = get_imFOV_mask_imgnd_xyz(cmmdl,[L_im,S_im],ddr_rov);

end
         

%%

DDRprj = [];
DDRprj.imxy = ddr_imxy;
DDRprj.NorthEastElevation = ddr_nee;
DDRprj.imFOV_mask = ddr_imFOV_mask;


end