function [msldemc_imFOVmask,msldemc_imFOVhdr,L_im,S_im,cmmdl_geo,option]...
    = mastcam_get_projMSLDEM2mastcam_v3_imFOVmask(MSLDEMdata,mastcamdata_obj,varargin)
% [msldemc_imFOVmask,msldemc_imFOVhdr,L_im,S_im,cmmdl_geo,option]...
%     = mastcam_get_projMSLDEM2mastcam_v3_imFOVmask(MSLDEMdata,mastcamdata_obj,varargin)
%   evaluate FOV of an image on an ortho-georeferenced image using a
%   georeferenced DEM image.
%  INPUTS:
%    MSLDEMdata: HSI class obj, 
%      MSLDEMdata at
%      https://astrogeology.usgs.gov/search/map/Mars/MarsScienceLaboratory/
%      Mosaics/MSL_Gale_DEM_Mosaic_10m
%    mastcamdata_obj: MASTCAMdata or MASTCAMgroup_eye class object
%  OUTPUTS
%   msldemc_imFOVmask: int8 image, [L_demc x S_demc]
%      2 if potentially in the FOV that have valid xy values
%      1 if potentially in the FOV but invalid xy values (in the wrong
%      direction)
%      0 if outside of the FOV
%      FOV should be true if the linear transformation shows the value
%      between -200 and 200+image_edge. The values outside of this are
%      considered to be inaccurate or invalid. In addition, pixels within
%      the distance 50m are considered to be in FOV, just in case.
%   msldemc_imFOVxy: [L_demc x S_demc x 2] 
%      page 1 is the x coordinate in the image frame,
%      page 2 is the y coordinate in the image frame.
%      x and y coordinate values outside of the range between -200 and 
%      200+image_edge are replaced with nans because the values outside of
%      this is far away from the calibration range, therefore accuracy is
%      not guaranteed. Size depend on FOV. The rectangle that minimally
%      encloses the FOV.
%   msldemc_imFOVhdr: struct having the fields below
%      msldemc_imFOVhdr.lines = len_vl;
%      msldemc_imFOVhdr.samples = len_vs;
%      msldemc_imFOVhdr.line_offset = line_offset;
%      msldemc_imFOVhdr.sample_offset = sample_offset;
%      msldemc_imFOVhdr.y = msldemc_northing;
%      msldemc_imFOVhdr.x = msldemc_easting;
%   L_im, S_im: scalar, size of the mastcam image 
%   cmmdl_geo: CAM_MDL class obj, 
%      camera model defined on the selected coordinate system  
%   option: struct, option of the method
%      field: COEF_MARGIN
%  OPTIONAL Parameters
%   "COORDINATE_SYSTEM": {'NEE','NORTHEASTNADIR','IAU_MARS_SPHERE',
%      'IAU_MARS_ELLIPSOID'}
%      coordinate system based on which the imFOVmask is calculated. 'NEE'
%      means north-east-elevation. (Actually, elevation is looking in the 
%      nadir direction in this coordinate system.)
%      (default) 'NORTHEASTNADIR'
%   "PROC_MODE": {'ENCLOSING_RECTANGLE',''SURROUNDING_COMPLEMENT'}
%      processing mode
%      (default) 'SURROUNDING_COMPLEMENT' 
%   "COEF_MARGIN"
%      (default) 2.1
%   


%% GET CAMERA, ROVER_NAV, and Image size information

coordsys = 'NorthEastNadir';
coef_mrgn = 2.1;
proc_mode = 'SURROUNDING_COMPLEMENT'; %'ENCLOSING_RECTANGLE';
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'COEF_MARGIN'
                coef_mrgn = varargin{i+1};
            case 'COORDINATE_SYSTEM'
                coordsys = varargin{i+1};
            case 'PROC_MODE'
                proc_mode = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end


cmmdl = mastcamdata_obj.CAM_MDL;
rover_nav_coord = mastcamdata_obj.ROVER_NAV;

switch upper(coordsys)
    case {'NEE','NORTHEASTELEVATION','NORTHEASTNADIR','NENADIR'}
        cmmdl_geo = transform_CAHVOR_MODEL_wROVER_NAV(cmmdl, ...
            rover_nav_coord);
        % cmmdl_geo.get_image_plane_unit_vectors();
    case {'IAU_MARS_SPHERE'}
        if ~isa(MSLDEMdata,'MSLGaleMosaicRadius_v3')
            error([ ...
                'MSLDEMdata needs to be an object of MSLGaleMosaicRadius_v3,' ...
                ' a subclass of MSLGaleDEMMosaic_v3' ...
                ]);
        end
        [cmmdl_geo] = transform_CAHVOR_MODEL_ROVERNAV2IAUMARS(cmmdl, ...
            rover_nav_coord,'Mars_Shape','Sphere');
        % Get mars radius
        [mars_radii] = mars_get_radii('Sphere');
        mars_re = mars_radii; mars_rp = mars_radii;
    case {'IAU_MARS_ELLIPSOID'}
        if ~isa(MSLDEMdata,'MSLGaleMosaicRadius_v3')
            error([ ...
                'MSLDEMdata needs to be an object of MSLGaleMosaicRadius_v3,' ...
                ' a subclass of MSLGaleDEMMosaic_v3' ...
                ]);
        end
        [cmmdl_geo] = transform_CAHVOR_MODEL_ROVERNAV2IAUMARS(cmmdl, ...
             rover_nav_coord,'Mars_Shape','Ellipsoid');
        [mars_radii] = mars_get_radii('Ellipsoid');
        mars_re = mars_radii(1); mars_rp = mars_radii(2);
    otherwise
        error('Undefined COORDINATE_SYSTEM %s',coordsys);
end
cmmdl_geo.get_image_plane_unit_vectors();

%-------------------------------------------------------------------------%
% Get the size of the mastcam image
%-------------------------------------------------------------------------%
L_im = mastcamdata_obj.L_im; S_im = mastcamdata_obj.S_im;

%% First compute dem_imFOV_mask
switch upper(coordsys)
    case {'NEE','NORTHEASTELEVATION','NORTHEASTNADIR'}
        switch mastcamdata_obj.Linearization
            case 1
                [msldem_imFOVmask] = mastcam_get_imFOVmask_msldem_cahvor_lr1(...
                    MSLDEMdata,L_im,S_im,cmmdl_geo,coef_mrgn,proc_mode);
            case 0
                [msldem_imFOVmask,lrange,srange]...
                    = mastcam_get_imFOVmask_msldem_cahvor_lr0(...
                    MSLDEMdata,L_im,S_im,cmmdl_geo,coef_mrgn,proc_mode);

            otherwise
                error('Linearization %d is not supported', ...
                    mastcamdata_obj.Linearization);
        end
    case {'IAU_MARS_SPHERE','IAU_MARS_ELLIPSOID'}
        switch mastcamdata_obj.Linearization
            case 1
                [msldem_imFOVmask] = mastcam_get_imFOVmask_msldem_cahvor_iaumars_lr1(...
                    MSLDEMdata,L_im,S_im,cmmdl_geo,coef_mrgn,proc_mode,mars_re,mars_rp);
            case 0
                [msldem_imFOVmask,lrange,srange]...
                    = mastcam_get_imFOVmask_msldem_cahvor_iaumars_lr0(...
                    MSLDEMdata,L_im,S_im,cmmdl_geo,coef_mrgn,proc_mode,mars_re,mars_rp);

            otherwise
                error('Linearization %d is not supported', ...
                    mastcamdata_obj.Linearization);
        end
        
    otherwise
        error('Undefined COORDINATE_SYSTEM %s',coordsys);
end

% [dem_imFOV_mask,dem_imFOV_mask_xyd] = get_imFOV_mask_MSLDEM(cmmdl_geo,...
%      [L_im,S_im],MSLDEMdata);

%% safeguarding
% the range of lines and samples are one pixel padded in all of the four
% directions.
msldem_imFOVmask_gt0 = msldem_imFOVmask>0;
valid_lines = find(any(msldem_imFOVmask_gt0',1));
lrnge = [max(valid_lines(1)-1,1), min(valid_lines(end)+1,MSLDEMdata.hdr.lines)];
len_vl = lrnge(2)-lrnge(1)+1;
valid_samples = find(any(msldem_imFOVmask_gt0,1));
srnge = [max(valid_samples(1)-1,1), min(valid_samples(end)+1,MSLDEMdata.hdr.samples)];
len_vs = srnge(2)-srnge(1)+1;

line_offset = lrnge(1)-1;
sample_offset = srnge(1)-1;

% cropping the image for the computation of sageguarding
msldemc_hdr_sg = [];
msldemc_hdr_sg.lines = len_vl;
msldemc_hdr_sg.samples = len_vs;
msldemc_hdr_sg.line_offset = line_offset;
msldemc_hdr_sg.sample_offset = sample_offset;



msldemc_imFOVmask = msldem_imFOVmask(lrnge(1):lrnge(2),srnge(1):srnge(2));
msldemc_imFOVmask(msldemc_imFOVmask<0) = 0;
clear msldem_imFOVmask;

switch upper(coordsys)
    case {'NEE','NORTHEASTELEVATION','NORTHEASTNADIR'}
        msldemc_northing = MSLDEMdata.northing(lrnge(1):lrnge(2));
        msldemc_easting  = MSLDEMdata.easting(srnge(1):srnge(2));
        tic; [safeguard_mask] = safeguard_msldemc_imFOVmask_napmc(...
            MSLDEMdata.imgpath,MSLDEMdata.hdr,msldemc_hdr_sg,...
            msldemc_northing,msldemc_easting,...
            msldemc_imFOVmask,cmmdl_geo); toc;
    case {'IAU_MARS_SPHERE','IAU_MARS_ELLIPSOID'}
        msldemc_latitude  = deg2rad(MSLDEMdata.latitude(lrnge(1):lrnge(2)));
        msldemc_longitude = deg2rad(MSLDEMdata.longitude(srnge(1):srnge(2)));
        tic; [safeguard_mask] = safeguard_iaumars_msldemc_imFOVmask_napmc(...
            MSLDEMdata.imgpath,MSLDEMdata.hdr,msldemc_hdr_sg, ...
            msldemc_latitude,msldemc_longitude,MSLDEMdata.OFFSET, ...
            msldemc_imFOVmask,cmmdl_geo); toc;
    otherwise
        error('Undefined COORDINATE_SYSTEM %s',coordsys);
end

if any(safeguard_mask,'all')
    fprintf('There seems to be a pixel that potentially in FOV and negatvie apmc');
end


msldemc_imFOVmask = msldemc_imFOVmask + safeguard_mask;

clear safeguard_mask; 


%% Next computing image coordinate for each pixel of DEM image.
% This computation is performed after the image mask is calculated because
% until evaluating the FOV, we do not know how to set valid_samples. We
% want to make this less memory intensive. 

valid_lines = find(any(msldemc_imFOVmask',1));
lrnge = [valid_lines(1),valid_lines(end)];
len_vl = lrnge(2)-lrnge(1)+1;
valid_samples = find(any(msldemc_imFOVmask,1));
srnge = [valid_samples(1),valid_samples(end)];
len_vs = srnge(2)-srnge(1)+1;

line_offset = lrnge(1)-1;
sample_offset = srnge(1)-1;

% msldemc_northing = msldemc_northing(lrnge(1):lrnge(2));
% msldemc_easting  = msldemc_easting(srnge(1):srnge(2));

msldemc_imFOVhdr = [];
msldemc_imFOVhdr.lines   = len_vl;
msldemc_imFOVhdr.samples = len_vs;
msldemc_imFOVhdr.line_offset   = line_offset + msldemc_hdr_sg.line_offset;
msldemc_imFOVhdr.sample_offset = sample_offset + msldemc_hdr_sg.sample_offset;
% msldemc_imFOVhdr.y = msldemc_northing;
% msldemc_imFOVhdr.x = msldemc_easting;

msldemc_imFOVmask = msldemc_imFOVmask(lrnge(1):lrnge(2),srnge(1):srnge(2));

%
switch mastcamdata_obj.Linearization
    case 1
        option = struct('COEF_MARGIN',coef_mrgn, ...
            'COORDINATE_SYSTEM',coordsys);
    case 0
        option = struct('COEF_MARGIN',coef_mrgn, ...
            'LRange',lrange,'SRange',srange,'PROC_MODE',proc_mode, ...
            'COORDINATE_SYSTEM',coordsys);
    otherwise
        error('Linearization %d is not supported',mastcamdata_obj.Linearization);
end


end