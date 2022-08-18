function [msldemc_imUFOVxynn,msldemc_imUFOVxy] = mastcam_calc_imUFOVxynn( ...
    mastcamdata_obj, MSLDEMdata, msldemc_imUFOVmask_obj,varargin)
% [msldemc_imUFOVxynn,msldemc_imUFOVxy] = mastcam_calc_imUFOVxynn( ...
%     mastcamdata_obj, MSLDEMdata, msldemc_imUFOVmask)
%  Get msldemc_imUFOVxynn (nearest neighbor xy pixels)
% INPUTS
%   mastcamdata_obj    : MASTCAMdata,MASTCAMgroup_wProcCode
%   MSLDEMdata         : MSLGaleDEMMosaic_v3
%   msldemc_imUFOVmask : ENVIRasterSingleLayerMSLDEMCProj
% OUTPUTS
%   msldemc_imUFOVxynn : 2 bands, the image of the same size as 
%   msldemc_imUFOVmask whose pixels are the rounded image xy coordinate 
%   values masked by msldemc_imUFOVmask.
%   msldemc_imUFOVxy : 2 bands, the image of the same size as
%   the msldemc_imUFOVmask whose pixels are the image xy coordinate values 
%   masked by msldemc_imUFOVmask.

coordsys = 'NorthEastNadir';

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})               
            % ## PROCESSING OPTIONS #--------------------------------------
            case {'COORDINATE_SYSTEM'}
                coordsys = varargin{i+1};
                validateattributes(coordsys,{'char'},{},mfilename,'COORDINATE_SYSTEM');
                
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

validateattributes(mastcamdata_obj,{'MASTCAMdata','MASTCAMgroup_wProcCode'}, ...
    {},mfilename,'mastcamdata_obj');
validateattributes(MSLDEMdata,{'MSLGaleDEMMosaic_v3'}, ...
    {},mfilename,'MSLDEMdata');
validateattributes(msldemc_imUFOVmask_obj, ...
    {'ENVIRasterMultBandMSLDEMCProj','ENVIRasterSingleLayerMSLDEMCProj'}, ...
    {},mfilename,'msldemc_imUFOVmask_obj');

%---------------------------------------------------------------------%
% Get cam_C_geo
%---------------------------------------------------------------------%
cmmdl = mastcamdata_obj.CAM_MDL;
rover_nav_coord = mastcamdata_obj.ROVER_NAV;
switch upper(coordsys)
    case {'NEE','NORTHEASTELEVATION','NORTHEASTNADIR','NENADIR'}
        cmmdl_geo = transform_CAHVOR_MODEL_wROVER_NAV(cmmdl, ...
            rover_nav_coord);
    case {'IAU_MARS_SPHERE'}
        if ~isa(MSLDEMdata,'MSLGaleMosaicRadius_v3')
            error([ ...
                'MSLDEMdata needs to be an object of MSLGaleMosaicRadius_v3,' ...
                ' a subclass of MSLGaleDEMMosaic_v3 for IAU_MARS_SPHERE' ...
                ]);
        end
        [cmmdl_geo] = transform_CAHVOR_MODEL_ROVERNAV2IAUMARS(cmmdl, ...
             rover_nav_coord,'Mars_Shape','Sphere');
    case {'IAU_MARS_ELLIPSOID'}
        if ~isa(MSLDEMdata,'MSLGaleMosaicRadius_v3')
            error([ ...
                'MSLDEMdata needs to be an object of MSLGaleMosaicRadius_v3,' ...
                ' a subclass of MSLGaleDEMMosaic_v3 for IAU_MARS_ELLIPSOID' ...
                ]);
        end
        [cmmdl_geo] = transform_CAHVOR_MODEL_ROVERNAV2IAUMARS(cmmdl, ...
             rover_nav_coord,'Mars_Shape','Ellipsoid');
    otherwise
        error('Undefined COORDINATE_SYSTEM %s',coordsys);
end
cmmdl_geo.get_image_plane_unit_vectors();

%%
if isempty(msldemc_imUFOVmask_obj.img)
    msldemc_imUFOVmask_img = msldemc_imUFOVmask_obj.readimg('precision','raw');
else
    msldemc_imUFOVmask_img = msldemc_imUFOVmask_obj.img;
end
srange = msldemc_imUFOVmask_obj.get_xrange_base();
lrange = msldemc_imUFOVmask_obj.get_yrange_base();

%%
switch upper(coordsys)
    case {'NEE','NORTHEASTELEVATION','NORTHEASTNADIR','NENADIR'}
        msldemc_northing = MSLDEMdata.northing(lrange(1):lrange(2));
        msldemc_easting  = MSLDEMdata.easting(srange(1):srange(2));
        switch mastcamdata_obj.Linearization
            case 1
                tic; [msldemc_imUFOVx,msldemc_imUFOVy] = cahv_get_imxy_MSLDEM_mex(...
                    MSLDEMdata.imgpath,MSLDEMdata.hdr,msldemc_imUFOVmask_obj.chdr,...
                    msldemc_northing,msldemc_easting,msldemc_imUFOVmask_img,cmmdl_geo); toc;
                % msldemc_imUFOVxy = cat(3,dem_imx,dem_imy);
            case 0
                tic; [msldemc_imUFOVx,msldemc_imUFOVy] = cahvor_get_imxy_msldem_mex(...
                    MSLDEMdata.imgpath,MSLDEMdata.hdr,msldemc_imUFOVmask_obj.chdr,...
                    msldemc_northing,msldemc_easting,msldemc_imUFOVmask_img,cmmdl_geo); toc;
                % msldemc_imUFOVxy = cat(3,dem_imx,dem_imy);

            otherwise
                error('Linearization %d is not supported',mastcamdata_obj.Linearization);
        end
    case {'IAU_MARS_SPHERE','IAU_MARS_ELLIPSOID'}
        msldemc_latitude  = deg2rad(MSLDEMdata.latitude(lrange(1):lrange(2)));
        msldemc_longitude = deg2rad(MSLDEMdata.longitude(srange(1):srange(2)));
        switch mastcamdata_obj.Linearization
            case 1
                tic; [msldemc_imUFOVx,msldemc_imUFOVy] = cahv_iaumars_get_imxy_MSLDEM_mex(...
                    MSLDEMdata.imgpath,MSLDEMdata.hdr, MSLDEMdata.OFFSET, ...
                    msldemc_imUFOVmask_obj.chdr,...
                    msldemc_latitude,msldemc_longitude,msldemc_imUFOVmask_img,cmmdl_geo); toc;
            case 0
                tic; [msldemc_imUFOVx,msldemc_imUFOVy] = cahvor_iaumars_get_imxy_MSLDEM_mex(...
                    MSLDEMdata.imgpath,MSLDEMdata.hdr, MSLDEMdata.OFFSET, ...
                    msldemc_imUFOVmask_obj.chdr,...
                    msldemc_latitude,msldemc_longitude,msldemc_imUFOVmask_img,cmmdl_geo); toc;

            otherwise
                error('Linearization %d is not supported',mastcamdata_obj.Linearization);
        end

    otherwise
        error('Undefined COORDINATE_SYSTEM %s',coordsys);
end

%% Get nearest neighbor map
% Mapping from msldemc to the nearest pixel of the MASTCAM image
msldemc_imUFOVxnn = round(msldemc_imUFOVx);
msldemc_imUFOVynn = round(msldemc_imUFOVy);

L_im = mastcamdata_obj.L_im; S_im = mastcamdata_obj.S_im;
msldemc_imUFOVxnn(msldemc_imUFOVxnn<0) = 0;
msldemc_imUFOVxnn(msldemc_imUFOVxnn>(S_im-1)) = S_im-1;
msldemc_imUFOVxnn(isnan(msldemc_imUFOVxnn)) = -1;
msldemc_imUFOVxnn = int16(msldemc_imUFOVxnn);
msldemc_imUFOVynn(msldemc_imUFOVynn<0) = 0;
msldemc_imUFOVynn(msldemc_imUFOVynn>(L_im-1)) = L_im-1;
msldemc_imUFOVynn(isnan(msldemc_imUFOVynn)) = -1;
msldemc_imUFOVynn = int16(msldemc_imUFOVynn);

msldemc_imUFOVxynn = cat(3,msldemc_imUFOVxnn,msldemc_imUFOVynn);

if nargout>1
    msldemc_imUFOVxy = cat(3,msldemc_imUFOVx,msldemc_imUFOVy);
end

end