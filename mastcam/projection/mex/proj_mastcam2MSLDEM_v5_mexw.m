function [mastcam_NEE,mastcam_msldemc_ref,mastcam_range,mastcam_msldemc_nn,...
    msldemc_imFOVmask_ctrnn,mastcam_emi,mastcam_surfplnc] = ...
    proj_mastcam2MSLDEM_v5_mexw(mastcamdata_obj,MSLDEMdata,MSTprj,varargin)
% proj_mastcam2MSLDEM_v5_mexw(mastcamdata_obj,MSLDEMdata,MSTprj)
%   Project mastcam image onto MSLDEMdata
%  INPUTS:
%    mastcamdata_obj: MASTCAMdata or MASTCAMgroup_eye class object
%    MSLDEMdata: HSI class obj, 
%      MSLDEMdata at
%      https://astrogeology.usgs.gov/search/map/Mars/MarsScienceLaboratory/
%      Mosaics/MSL_Gale_DEM_Mosaic_10m
%    MSTprj: object of class MASTCAMCameraProjectionMSLDEM
%  OUTPUTS:
%     mastcam_NEE: [L_im x S_im x 3]
%       pages 1,2,3 are northing, easting, and elevation. 
%     mastcam_msldemc_ref: [L_im x S_im x 3]
%       indicate which triangle in the MSLDEMdata, the pixel is located.
%       pages 1,2,3 are sample, line, and j. j is either 1 or 2, indicating
%       which triangle at the (sample, line).
%     mastcam_range: [L_im x S_im]
%       range for each pixel.
%     mastcam_msldemc_nn: [L_im x S_im x 2]
%       nearest neighbor pixels in DDR image. The first page is the sample
%       indices and the second is the line indexes.
%     msldemc_imFOVmask_ctrnn: [L_demc x S_demc x 1]
%       Boolean, imFOV_mask with hidden points are removed.
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
l1 = MSTprj.msldemc_imFOVhdr.line_offset+1;
lend = MSTprj.msldemc_imFOVhdr.line_offset+MSTprj.msldemc_imFOVhdr.lines;
s1 = MSTprj.msldemc_imFOVhdr.sample_offset+1;
send = MSTprj.msldemc_imFOVhdr.sample_offset+MSTprj.msldemc_imFOVhdr.samples;


% dem_imFOVd_mask_crop = ~isnan(MSTprj.msldemc_imxy(:,:,1));

%% Main computation.
switch upper(coordsys)
    case {'NEE','NORTHEASTELEVATION','NORTHEASTNADIR','NENADIR'}
        dem_northing_crop = MSLDEMdata.northing(l1:lend);
        dem_easting_crop  = MSLDEMdata.easting(s1:send);
        switch mastcamdata_obj.Linearization
            case 1
                tic; [im_north,im_east,im_elev,msldem_refx,msldem_refy, ...
                    msldem_refs,im_range,im_nnx,im_nny,im_emi, ...
                    im_pnx,im_pny,im_pnz,im_pc] = ...
                    cahv_proj_mastcam2MSLDEM_v5_mex(...
                        MSLDEMdata.imgpath,...0
                        MSLDEMdata.hdr,...1
                        MSTprj.msldemc_imFOVhdr,...2
                        dem_northing_crop,...3
                        dem_easting_crop,...4
                        MSTprj.msldemc_imFOVmask,...5
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

                tic; [im_north,im_east,im_elev,msldem_refx,msldem_refy, ...
                    msldem_refs,im_range,im_nnx,im_nny,im_emi, ...
                    im_pnx,im_pny,im_pnz,im_pc] = ...
                    cahvor_proj_mastcam2MSLDEM_v5_mex(...
                        MSLDEMdata.imgpath,...0
                        MSLDEMdata.hdr,...1
                        MSTprj.msldemc_imFOVhdr,...2
                        dem_northing_crop,...3
                        dem_easting_crop,...4
                        MSTprj.msldemc_imFOVmask,...5
                        S_im,L_im,...6,7
                        cmmdl_geo,...8
                        PmC(:,:,1),PmC(:,:,2),PmC(:,:,3)); toc; % 9, 10, 11
                mastcam_range = sqrt(im_range);

            otherwise
                error('Linearization %d is not supported',mastcamdata_obj.Linearization);
        end
    case {'IAU_MARS_SPHERE','IAU_MARS_ELLIPSOID'}
        demc_latitude  = deg2rad(MSLDEMdata.latitude(l1:lend));
        demc_longitude = deg2rad(MSLDEMdata.longitude(s1:send));
        switch mastcamdata_obj.Linearization
            case 1
                tic; [im_north,im_east,im_elev,msldem_refx,msldem_refy, ...
                    msldem_refs,im_range,im_nnx,im_nny,im_emi, ...
                    im_pnx,im_pny,im_pnz,im_pc] = ...
                    cahv_iaumars_proj_mastcam2MSLDEM_v5_mex(...
                        MSLDEMdata.imgpath,         ...0
                        MSLDEMdata.hdr,             ...1
                        MSTprj.msldemc_imFOVhdr,    ...2
                        demc_latitude,              ...3
                        demc_longitude,             ...4
                        MSLDEMdata.OFFSET,          ...5
                        MSTprj.msldemc_imFOVmask,   ...6
                        S_im,L_im,                  ...7,8
                        cmmdl_geo,                  ...9
                        PmC(:,:,1),PmC(:,:,2),PmC(:,:,3)); toc; % 10,11,12
                mastcam_range = sqrt(im_range);
            case 0
                tic; [im_north,im_east,im_elev,msldem_refx,msldem_refy, ...
                    msldem_refs,mastcam_range,im_nnx,im_nny,im_emi, ...
                    im_pnx,im_pny,im_pnz,im_pc] = ...
                    cahvor_iaumars_proj_mastcam2MSLDEM_v6_mex(...
                        MSLDEMdata.imgpath,         ...0
                        MSLDEMdata.hdr,             ...1
                        MSTprj.msldemc_imFOVhdr,    ...2
                        demc_latitude,              ...3
                        demc_longitude,             ...4
                        MSLDEMdata.OFFSET,          ...5
                        MSTprj.msldemc_imFOVmask,   ...6
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
% Evaluate the neaerest neighbor from the center of each image pixels.
% im_nnx and im_nny are the nearest neighbor indices 
% (the first index starts with 0)
idx_ctrnn = im_nnx*MSTprj.msldemc_imFOVhdr.lines + (im_nny+1);
idx_ctrnn = idx_ctrnn(:);
idx_ctrnn_1d_nisnan = (idx_ctrnn>=1); % invalid pixels are filled with -1.
idx_ctrnn_1d_ok = idx_ctrnn(idx_ctrnn_1d_nisnan);
img_mask_ctrnn = false(MSTprj.msldemc_imFOVhdr.lines*MSTprj.msldemc_imFOVhdr.samples,1);
img_mask_ctrnn(idx_ctrnn_1d_ok) = true;
img_mask_ctrnn = reshape(img_mask_ctrnn,[MSTprj.msldemc_imFOVhdr.lines,MSTprj.msldemc_imFOVhdr.samples]);
img_mask_ctrnn = int8(img_mask_ctrnn);

%% Fill propoerties of MASTCAMCameraProjectionMSLDEM object.
mastcam_NEE = cat(3,im_north,im_east,im_elev);
mastcam_msldemc_ref = cat(3,msldem_refx,msldem_refy,msldem_refs); % original, index start from 0.
% mastcam_range = sqrt(im_range);
mastcam_msldemc_nn = cat(3,im_nnx+1,im_nny+1); % index start from 1.
msldemc_imFOVmask_ctrnn = img_mask_ctrnn;
mastcam_emi = acosd(im_emi);
mastcam_surfplnc = cat(3,im_pnx,im_pny,im_pnz,im_pc);

end
