function [mastcam_XYZ,mastcam_latlon,mastcam_NE,mastcam_zenith, ...
    mastcam_ref_msldem,mastcam_range,mastcam_nn_msldem,...
    msldemc_imFOVmask_pxlctrnn,mastcam_emi,mastcam_surfplnc,option] = ...
    proj_mastcam2MSLDEM_v5_mexw_v2(mastcamdata_obj,MSLDEMdata,msldemc_imFOVmask,varargin)
% proj_mastcam2MSLDEM_v5_mexw_v2(mastcamdata_obj,MSLDEMdata,MSTprj)
%   Project mastcam image onto MSLDEMdata
%  INPUTS:
%    mastcamdata_obj: MASTCAMdata or MASTCAMgroup_eye class object
%    MSLDEMdata: MSLGaleDEMMosaic_v3 class obj, 
%      MSLDEMdata at
%      https://astrogeology.usgs.gov/search/map/Mars/MarsScienceLaboratory/
%      Mosaics/MSL_Gale_DEM_Mosaic_10m
%    msldemc_imFOVmask: object of class ENVIRasterSingleLayerMSLDEMCProj
%  OUTPUTS:
%     mastcam_XYZ: [L_im x S_im x 3] ,double
%       pages 1,2 are X, Y, and Z. (With NorthEastNadir, XYZ are northing,
%       easting, nadir (negative direction of the elevation. With IAU_MARS,
%       X, Y, and Z are XYZ of the rectangular coordinate system.
%     mastcam_latlon: [L_im x S_im x 2]
%       pages 1,2 are latitude and longitude in degree.
%     mastcam_NE: [L_im x S_im x 2]
%       pages 1,2 are northing, easting.
%     mastcam_zenith: [L_im x S_im x 1]
%       the coordinate value in the zenith direction. It depends on the
%       input. Radius with IAU_MARS, and Elevation with "North-East-Nadir".
%     mastcam_ref_msldemc: [L_im x S_im x 2]
%       indicate which triangle in the MSLDEMdata, the pixel is located.
%       pages 1,2,3 are sample, line, and j. j is either 0 or 1, indicating
%       which triangle at the (sample, line).
%     mastcam_range: [L_im x S_im]
%       range for each pixel.
%     mastcam_nn_msldem: [L_im x S_im x 2]
%       nearest neighbor pixels in MSLDEM image. The first page is the sample
%       indices and the second is the line indexes of 
%     msldemc_imFOVmask_pxlctrnn: [L_demc x S_demc x 1]
%       Boolean, imFOV_mask with hidden points are removed.
%     option: processing options
%     mastcam_emi: [L_im x S_im]
%       double, emission angles at each pixel center
%     mastcam_surfplnc: [L_im x S_im x 4]
%       double, parameters for the surface 
%       
%  OPTIONAL Parameters
%   "COORDINATE_SYSTEM": {'NEE','NORTHEASTNADIR','IAU_MARS_SPHERE',
%      'IAU_MARS_ELLIPSOID'}
%      coordinate system based on which the imFOVmask is calculated. 'NEE'
%      means north-east-elevation. (Actually, elevation is looking in the 
%      nadir direction in this coordinate system.)
%      (default) 'NORTHEASTNADIR'
%
% Copyright (C) 2021 Yuki Itoh <yukiitohand@gmail.com>

coordsys = 'NorthEastNadir';
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'COORDINATE_SYSTEM'
                coordsys = varargin{i+1};
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
    case {'IAU_MARS_ELLIPSOID'}
        if ~isa(MSLDEMdata,'MSLGaleMosaicRadius_v3')
            error([ ...
                'MSLDEMdata needs to be an object of MSLGaleMosaicRadius_v3,' ...
                ' a subclass of MSLGaleDEMMosaic_v3' ...
                ]);
        end
        [cmmdl_geo] = transform_CAHVOR_MODEL_ROVERNAV2IAUMARS(cmmdl, ...
             rover_nav_coord,'Mars_Shape','Ellipsoid');
    otherwise
        error('Undefined COORDINATE_SYSTEM %s',coordsys);
end
cmmdl_geo.get_image_plane_unit_vectors();


%-------------------------------------------------------------------------%
% Get the size of the mastcam image
%-------------------------------------------------------------------------%
L_im = mastcamdata_obj.L_im; S_im = mastcamdata_obj.S_im;

%-------------------------------------------------------------------------%
% Construct Image grid
%-------------------------------------------------------------------------%
imx_im_1d = 0:(S_im-1);
imy_im_1d = reshape(0:(L_im-1),[],1);
[imx_im,imy_im] = meshgrid(imx_im_1d,imy_im_1d);
im_imxy_vec2d = permute(reshape(cat(2,imx_im,imy_im),[L_im*S_im,2]),[2,1,3]);

PmC = cmmdl_geo.get_p_minus_c_from_xy(im_imxy_vec2d);
PmC = reshape(PmC',[L_im,S_im,3]);


%%
%-------------------------------------------------------------------------%
% Get MSL DEM information
%-------------------------------------------------------------------------%
if isempty(msldemc_imFOVmask.img)
    msldemc_imFOVmask_img = msldemc_imFOVmask.readimg('precision','raw');
else
    msldemc_imFOVmask_img = msldemc_imFOVmask.img;
end
srange = msldemc_imFOVmask.get_xrange_base();
lrange = msldemc_imFOVmask.get_yrange_base();
% l1   = msldemc_imFOVmask.get_y_base(1);
% lend = msldemc_imFOVmask.get_y_base(msldemc_imFOVmask.hdr.lines);
% s1   = msldemc_imFOVmask.get_x_base(1);
% send = msldemc_imFOVmask.get_x_base(msldemc_imFOVmask.hdr.samples);


% dem_imFOVd_mask_crop = ~isnan(MSTprj.msldemc_imxy(:,:,1));

%% Main computation.
switch upper(coordsys)
    case {'NEE','NORTHEASTELEVATION','NORTHEASTNADIR','NENADIR'}
        dem_northing_crop = MSLDEMdata.northing(lrange(1):lrange(2));
        dem_easting_crop  = MSLDEMdata.easting(srange(1):srange(2));
        switch mastcamdata_obj.Linearization
            case 1
                tic; [im_x,im_y,im_z,msldemc_refx,msldemc_refy, ...
                    msldemc_refs,im_range,im_nnx,im_nny,im_emi, ...
                    im_pnx,im_pny,im_pnz,im_pc] = ...
                    cahv_proj_mastcam2MSLDEM_v5_mex(...
                        MSLDEMdata.imgpath,...0
                        MSLDEMdata.hdr,...1
                        msldemc_imFOVmask.chdr,...2
                        dem_northing_crop,...3
                        dem_easting_crop,...4
                        msldemc_imFOVmask_img,...5
                        S_im,L_im,...6,7
                        cmmdl_geo,...8
                        PmC(:,:,1),PmC(:,:,2),PmC(:,:,3)); toc; % 9, 10, 11
                mastcam_range = sqrt(im_range);
                % comment out below is for an obsolete function. superseded by the
                % function above because it is faster and memory efficient.
                % tic; [im_north,im_east,im_elev,msldem_refx,msldem_refy,msldem_refs,im_range,...
                %         im_nnx,im_nny,im_emi,im_pnx,im_pny,im_pnz,im_pc] = ...
                %     proj_mastcam2MSLDEM_v4_mex(...
                %         MSLDEMdata.imgpath,...0
                %         MSLDEMdata.hdr,...1
                %         MSTprj.msldemc_imFOVhdr,...2
                %         dem_northing_crop,...3
                %         dem_easting_crop,...4
                %         MSTprj.msldemc_imFOVxy(:,:,1),...5
                %         MSTprj.msldemc_imFOVxy(:,:,2),...6
                %         MSTprj.msldemc_imFOVmask,...7
                %         S_im,L_im,...8,9
                %         cmmdl_geo,...10
                %         PmC(:,:,1),PmC(:,:,2),PmC(:,:,3)); toc; % 11,12,13

            case 0

                tic; [im_x,im_y,im_z,msldemc_refx,msldemc_refy, ...
                    msldemc_refs,im_range,im_nnx,im_nny,im_emi, ...
                    im_pnx,im_pny,im_pnz,im_pc] = ...
                    cahvor_proj_mastcam2MSLDEM_v5_mex(...
                        MSLDEMdata.imgpath,...0
                        MSLDEMdata.hdr,...1
                        msldemc_imFOVmask.chdr,...2
                        dem_northing_crop,...3
                        dem_easting_crop,...4
                        msldemc_imFOVmask_img,...5
                        S_im,L_im,...6,7
                        cmmdl_geo,...8
                        PmC(:,:,1),PmC(:,:,2),PmC(:,:,3)); toc; % 9, 10, 11
                mastcam_range = sqrt(im_range);

            otherwise
                error('Linearization %d is not supported',mastcamdata_obj.Linearization);
        end
    case {'IAU_MARS_SPHERE','IAU_MARS_ELLIPSOID'}
        demc_latitude  = deg2rad(MSLDEMdata.latitude(lrange(1):lrange(2)));
        demc_longitude = deg2rad(MSLDEMdata.longitude(srange(1):srange(2)));
        switch mastcamdata_obj.Linearization
            case 1
                tic; [im_x,im_y,im_z,msldemc_refx,msldemc_refy, ...
                    msldemc_refs,im_range,im_nnx,im_nny,im_emi, ...
                    im_pnx,im_pny,im_pnz,im_pc] = ...
                    cahv_iaumars_proj_mastcam2MSLDEM_v5_mex(...
                        MSLDEMdata.imgpath,         ...0
                        MSLDEMdata.hdr,             ...1
                        msldemc_imFOVmask.chdr,    ...2
                        demc_latitude,              ...3
                        demc_longitude,             ...4
                        MSLDEMdata.OFFSET,          ...5
                        msldemc_imFOVmask_img,   ...6
                        S_im,L_im,                  ...7,8
                        cmmdl_geo,                  ...9
                        PmC(:,:,1),PmC(:,:,2),PmC(:,:,3)); toc; % 10,11,12
                mastcam_range = sqrt(im_range);
            case 0
                tic; [im_x,im_y,im_z,msldemc_refx,msldemc_refy, ...
                    msldemc_refs,mastcam_range,im_nnx,im_nny,im_emi, ...
                    im_pnx,im_pny,im_pnz,im_pc] = ...
                    cahvor_iaumars_proj_mastcam2MSLDEM_v6_mex(...
                        MSLDEMdata.imgpath,         ...0
                        MSLDEMdata.hdr,             ...1
                        msldemc_imFOVmask.chdr,    ...2
                        demc_latitude,              ...3
                        demc_longitude,             ...4
                        MSLDEMdata.OFFSET,          ...5
                        msldemc_imFOVmask_img,   ...6
                        S_im,L_im,                  ...7,8
                        cmmdl_geo,                  ...9
                        PmC(:,:,1),PmC(:,:,2),PmC(:,:,3)); toc; % 10,11,12
                
            otherwise
                error('Linearization %d is not supported',mastcamdata_obj.Linearization);
        end  
        
    otherwise
        error('Undefined COORDINATE_SYSTEM %s',coordsys);
end



%% Post computation task.

% First create
% Evaluate the neaerest neighbor from the center of each image pixels.
% im_nnx and im_nny are the nearest neighbor indices 
% (the first index starts with 0)
idx_ctrnn = im_nnx*msldemc_imFOVmask.hdr.lines + (im_nny+1);
idx_ctrnn = idx_ctrnn(:);
idx_ctrnn_1d_nisnan = (idx_ctrnn>=1); % invalid pixels are filled with -1.
idx_ctrnn_1d_ok = idx_ctrnn(idx_ctrnn_1d_nisnan);
msldemc_imFOVmask_pxlctrnn = false(msldemc_imFOVmask.hdr.lines*msldemc_imFOVmask.hdr.samples,1);
msldemc_imFOVmask_pxlctrnn(idx_ctrnn_1d_ok) = true;
msldemc_imFOVmask_pxlctrnn = reshape(msldemc_imFOVmask_pxlctrnn,[msldemc_imFOVmask.hdr.lines,msldemc_imFOVmask.hdr.samples]);
msldemc_imFOVmask_pxlctrnn = int8(msldemc_imFOVmask_pxlctrnn);

%
sample_offset = msldemc_imFOVmask.chdr.sample_offset;
line_offset   = msldemc_imFOVmask.chdr.line_offset;
mastcam_XYZ = cat(3,im_x,im_y,im_z);

% reference triangle, compensate offset and convert the index style to MATLAB
msldem_refx = msldemc_refx; msldem_refy = msldemc_refy; msldem_refs = msldemc_refs;
ref_valid = msldem_refx>-1;
msldem_refx(ref_valid) = msldem_refx(ref_valid)+(1+sample_offset); % MATLAB style, compensate offset
msldem_refy(ref_valid) = msldem_refy(ref_valid)+(1+line_offset);   % MATLAB style, compensate offset
msldem_refs(ref_valid) = msldem_refs(ref_valid) + 1;               % MATLAB style
mastcam_ref_msldem = cat(3,msldem_refx,msldem_refy,msldem_refs); 
%index start from 1.

% nearest neighbor, compensate offset and convert the index style to MATLAB
nnx_msldem = im_nnx; nny_msldem = im_nny; ref_valid = nnx_msldem>-1;
nnx_msldem(ref_valid) = nnx_msldem(ref_valid) + (1+sample_offset);
nny_msldem(ref_valid) = nny_msldem(ref_valid) + (1+line_offset);
mastcam_nn_msldem = cat(3,nnx_msldem,nny_msldem);
% index start from 1.

mastcam_emi = acosd(im_emi);
mastcam_surfplnc = cat(3,im_pnx,im_pny,im_pnz,im_pc);

switch upper(coordsys)
    case {'NEE','NORTHEASTELEVATION','NORTHEASTNADIR','NENADIR'}
        mastcam_NE = mastcam_XYZ(:,:,1:2);
        mastcam_lat = MSLDEMdata.latitude(MSLDEMdata.northing2y(mastcam_XYZ(:,:,1)));
        mastcam_lon = MSLDEMdata.longitude(MSLDEMdata.easting2x(mastcam_XYZ(:,:,2)));
        mastcam_latlon = cat(3,mastcam_lat,mastcam_lon);
        mastcam_zenith = mastcam_XYZ(:,:,3); % zenith is the elevation in this case.
    case {'IAU_MARS_SPHERE','IAU_MARS_ELLIPSOID'}
        mastcam_zenith = sqrt(sum(mastcam_XYZ.^2,3)); % zenith is the radius in this case.
        mastcam_lat = asind(mastcam_XYZ(:,:,3)./mastcam_zenith);
        mastcam_lon = atan2d(mastcam_XYZ(:,:,2), mastcam_XYZ(:,:,1));
        mastcam_latlon   = cat(3,mastcam_lat,mastcam_lon);
        mastcam_northing = MSLDEMdata.northing(MSLDEMdata.lat2y(mastcam_lat));
        mastcam_easting  = MSLDEMdata.easting(MSLDEMdata.lon2x(mastcam_lon));
        mastcam_NE = cat(3,mastcam_northing,mastcam_easting);
    otherwise
        error('Undefined COORDINATE_SYSTEM %s',coordsys);
end


option = [];
switch upper(coordsys)
    case {'NEE','NORTHEASTELEVATION','NORTHEASTNADIR','NENADIR'}
        option.xyz_coordsys = 'NorthEastNadir';
    case {'IAU_MARS_SPHERE','IAU_MARS_ELLIPSOID'}
        option.xyz_coordsys = 'IAU_MARS_SPHERE';
    otherwise
        error('Undefined COORDINATE_SYSTEM %s',coordsys);
end

option.linearization = mastcamdata_obj.Linearization;
option.cmmdl_geo     = cmmdl_geo;
option.msldemc_imFOVmask = msldemc_imFOVmask.basename;



end
