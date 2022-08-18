function [msldemc_imUFOVmask] = mastcam_getUFOVwMSLDEM_mexw_v2( ...
    mastcamdata_obj,MSLDEMdata,msldemc_imFOVmask,msldemc_imFOVmask_pxlctrnn_obj,varargin)
%[msldemc_imUFOVmask] = mastcam_getUFOVwMSLDEM_mexw(MSTproj,varargin)
%   calculate image Unobstructed Field of View (UFOV) with the precomputed 
%   imageFOV mask. 
%  INPUTS
%   mastcamdata_obj: either of 'MASTCAMdata' or 'MASTCAMgroup_wProcCode'
%   MSLDEMdata : MSLDEMGaleMosaic_v3
%   msldemc_imFOVmask_obj : ENVIRasterSingleLayerMSLDEMCProj
%   msldemc_imFOVmask_pxlctrnn_obj: ENVIRasterSingleLayerMSLDEMCProj
%  OUTPUTS
%   msldemc_imUFOVmask: int8 2 dimensonal array. The size is defined in 
%   msldemc_imFOVmask
%  OPTIONAL PARAMETERS
%   "UPDATE_CACHE": boolean, force updating cache files or not
%     (default) false
%    'CACHE_DIRPATH': 
%       (default) msl_env_vars.dirpath_cache
%    'SAVE_FILE': 
%       (default) 1
%    'FORCE': boolean, whether to process forcefully
%       (default) 0
%    {'LOAD_CACHE_IFEXIST','LCIE'}
%       (default) 1
%   "BORDER_ASSESS_OPT": string, combination of 'dcl'
%     (default) 'd'
%   "COORDINATE_SYSTEM": {'NorthEastNadir','IAU_MARS_SPHERE',IAU_MARS_ELLIPSOID'}
%     (default) 'NorthEastNadir'
%   "K_L","K_S": double, scaling coefficients in the line and sample
%     directions for rectangle bins
%     (default) 1 (no scaling, bin size is 1 x 1 [pxl])
%   "PROC_MODE": char array, processing mode. See other documents for the
%   meaning.
%     (default) 'L2PBK_LL0DYU_M3'

global msl_env_vars
pdir_cache = msl_env_vars.dirpath_cache;
coordsys = 'NorthEastNadir';
border_assess_opt = 'd';
save_file = true;
force = false;
load_cache_ifexist = true;
cache_vr = ''; % {'v0','v1'}
PROC_MODE = 'L2PBK_LL0DYU_M3';
K_L = 1; K_S = 1;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            % ## I/O OPTIONS #---------------------------------------------
            case 'CACHE_DIRPATH'
                pdir_cache = varargin{i+1};
                validateattributes(pdir_cache,{'char'},{},mfilename,'CACHE_DIRPATH');
            case 'SAVE_FILE'
                save_file = varargin{i+1};
                validateattributes(save_file,{'numeric','logical'},{'binary'},mfilename,'SAVE_FILE');
            case 'FORCE'
                force = varargin{i+1};
                validateattributes(force,{'numeric','logical'},{'binary'},mfilename,'FORCE');
            case {'LOAD_CACHE_IFEXIST','LCIE'}
                load_cache_ifexist = varargin{i+1};
                validateattributes(load_cache_ifexist,{'numeric','logical'}, ...
                    {'binary'},mfilename,'LOAD_CACHE_IFEXIST');
            case {'CACHE_VER','CACHE_VERSION'}
                cache_vr = varargin{i+1};
                validateattributes(cache_vr,{'char'},{},mfilename,'CACHE_VERSION');
            
            % ## PROCESSING OPTIONS #--------------------------------------
            case 'BORDER_ASSESS_OPT'
                border_assess_opt = lower(varargin{i+1});
                
            case 'COORDINATE_SYSTEM'
                coordsys = varargin{i+1};
                
            case 'PROC_MODE'
                PROC_MODE = varargin{i+1};
                
                
            case 'K_L'
                K_L = varargin{i+1};
            case 'K_S'
                K_S = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

if isempty(cache_vr)
   error('Please enter cache version with the option "CACHE_VER" or "CACHE_VERSION".'); 
end


% mastcamdata_obj = MSTproj.MASTCAMdata;
% MSLDEMdata = MSTproj.MSLDEMdata;
%-------------------------------------------------------------------------%
% Get the size of the mastcam image
%-------------------------------------------------------------------------%
L_im = mastcamdata_obj.L_im; S_im = mastcamdata_obj.S_im;

%-------------------------------------------------------------------------%
% Get cam_C_geo
%-------------------------------------------------------------------------%
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
                ' a subclass of MSLGaleDEMMosaic_v3 for IAU_MARS_SPHERE'      ...
                ]);
        end
        [cmmdl_geo] = transform_CAHVOR_MODEL_ROVERNAV2IAUMARS(cmmdl, ...
             rover_nav_coord,'Mars_Shape','Sphere');
    case {'IAU_MARS_ELLIPSOID'}
        if ~isa(MSLDEMdata,'MSLGaleMosaicRadius_v3')
            error([ ...
                'MSLDEMdata needs to be an object of MSLGaleMosaicRadius_v3,' ...
                ' a subclass of MSLGaleDEMMosaic_v3 for IAU_MARS_ELLIPSOID'   ...
                ]);
        end
        [cmmdl_geo] = transform_CAHVOR_MODEL_ROVERNAV2IAUMARS(cmmdl, ...
             rover_nav_coord,'Mars_Shape','Ellipsoid');
    otherwise
        error('Undefined COORDINATE_SYSTEM %s',coordsys);
end
cmmdl_geo.get_image_plane_unit_vectors();
C_geo = cmmdl_geo.C;

%-------------------------------------------------------------------------%
% Create label name for cache files
%-------------------------------------------------------------------------%
[cache_dirname] = mastcam_get_cache_dirname(mastcamdata_obj);
dirpath_cache = joinPath(pdir_cache,cache_dirname);
if ~exist(dirpath_cache,'dir')
    mkdir(dirpath_cache);
end
[basename_cache_com] = mastcam_create_basename_cache(mastcamdata_obj);


%-------------------------------------------------------------------------%
% msldemc_imFOVmask
%-------------------------------------------------------------------------%
srange = msldemc_imFOVmask.get_xrange_base();
lrange = msldemc_imFOVmask.get_yrange_base();

if isempty(msldemc_imFOVmask.img)
    msldemc_imFOVmask_img = msldemc_imFOVmask.readimg('precision','raw');
else
    msldemc_imFOVmask_img = msldemc_imFOVmask.img;
end

%-------------------------------------------------------------------------%
% msldemc_imFOVmask_pxlctrnn_obj
%-------------------------------------------------------------------------%
if isempty(msldemc_imFOVmask_pxlctrnn_obj.img)
    msldemc_imFOVmask_pxlctrnn_img = msldemc_imFOVmask_pxlctrnn_obj.readimg('precision','raw');
else
    msldemc_imFOVmask_pxlctrnn_img = msldemc_imFOVmask_pxlctrnn_obj.img;
end

%%
%-------------------------------------------------------------------------%
% Evaluate the center of each pixel
%-------------------------------------------------------------------------%
basename_imUFOVmask_ctr = mastcam_get_basename_cache(basename_cache_com, ...
    'imUFOVmask_ctr',cache_vr);
[ufov_ctr_imgpath,ufov_ctr_hdrpath] = mastcam_get_cachefilepath( ...
    basename_imUFOVmask_ctr,dirpath_cache);

[flg_reproc] = doyouwanttoprocess({ufov_ctr_imgpath,ufov_ctr_hdrpath},force,load_cache_ifexist);
if flg_reproc
    
    switch upper(coordsys)
        case {'NEE','NORTHEASTELEVATION','NORTHEASTNADIR','NENADIR'}
            msldemc_northing = MSLDEMdata.northing(lrange(1):lrange(2));
            msldemc_easting  = MSLDEMdata.easting(srange(1):srange(2));
            % msldemc_xmc      = msldemc_northing - C_geo(1);
            % msldemc_ymc      = msldemc_easting  - C_geo(2);
            % tic; [ msldemc_imUFOVmask_ctr ] = get_msldemtUFOVmask_ctr_wmsldemc_L2_mex(...
            %     MSLDEMdata.imgpath,MSLDEMdata.hdr,msldemc_imFOVmask.chdr,...
            %     msldemc_xmc,msldemc_ymc,msldemc_imFOVmask_img,S_im,L_im,cmmdl_geo); toc;
            switch upper(PROC_MODE)
                case 'L2PBK_LL0DYU_M3'
                     tic; [ msldemc_imUFOVmask_ctr] =  ...
                        northeast_get_msldemUFOVmask_ctr_L2PBK_LL0_M3_mex(...
                        MSLDEMdata.imgpath,MSLDEMdata.hdr,  ...
                        msldemc_imFOVmask.chdr,...
                        msldemc_northing,msldemc_easting, ...
                        msldemc_imFOVmask_img,S_im,L_im,cmmdl_geo,K_L,K_S,1); toc;
                otherwise
                    error('Undefined PROC_MODE %s for %s',PROC_MODE,coordsys);
            end
            
        case {'IAU_MARS_SPHERE','IAU_MARS_ELLIPSOID'}
            msldemc_latitude  = deg2rad(MSLDEMdata.latitude(lrange(1):lrange(2)));
            msldemc_longitude = deg2rad(MSLDEMdata.longitude(srange(1):srange(2)));
            switch upper(PROC_MODE)
                case 'L2PBN_DAP_M2'
                    tic; [ msldemc_imUFOVmask_ctr] =  ...
                        iaumars_get_msldemtUFOVmask_ctr_L2PBN_DAP_M2_mex(...
                        MSLDEMdata.imgpath,MSLDEMdata.hdr, MSLDEMdata.OFFSET, ...
                        msldemc_imFOVmask.chdr,...
                        msldemc_latitude,msldemc_longitude, ...
                        msldemc_imFOVmask_img,S_im,L_im,cmmdl_geo); toc;
                case 'L2PBK_DAPDYM_M2'
                    tic; [ msldemc_imUFOVmask_ctr] =  ...
                        iaumars_get_msldemtUFOVmask_ctr_L2PBK_DAP_M2_mex(...
                        MSLDEMdata.imgpath,MSLDEMdata.hdr, MSLDEMdata.OFFSET, ...
                        msldemc_imFOVmask.chdr,...
                        msldemc_latitude,msldemc_longitude, ...
                        msldemc_imFOVmask_img,S_im,L_im,cmmdl_geo,K_L,K_S); toc;
                case 'L2PBK_DAADYM_M3'
                    tic; [ msldemc_imUFOVmask_ctr] =  ...
                        iaumars_get_msldemtUFOVmask_ctr_L2PBK_DAA_M3_mex(...
                        MSLDEMdata.imgpath,MSLDEMdata.hdr, MSLDEMdata.OFFSET, ...
                        msldemc_imFOVmask.chdr,...
                        msldemc_latitude,msldemc_longitude, ...
                        msldemc_imFOVmask_img,S_im,L_im,cmmdl_geo,K_L,K_S); toc;
                    
                case 'L2PBK_DAR_M3'
                    tic; [ msldemc_imUFOVmask_ctr] =  ...
                        iaumars_get_msldemtUFOVmask_ctr_L2PBK_DAR_M3_mex(...
                        MSLDEMdata.imgpath,MSLDEMdata.hdr, MSLDEMdata.OFFSET, ...
                        msldemc_imFOVmask.chdr,...
                        msldemc_latitude,msldemc_longitude, ...
                        msldemc_imFOVmask_img,S_im,L_im,cmmdl_geo,K_L,K_S,0); toc;
                case 'L2PBK_DARDYM_M3'
                    tic; [ msldemc_imUFOVmask_ctr] =  ...
                        iaumars_get_msldemtUFOVmask_ctr_L2PBK_DAR_M3_mex(...
                        MSLDEMdata.imgpath,MSLDEMdata.hdr, MSLDEMdata.OFFSET, ...
                        msldemc_imFOVmask.chdr,...
                        msldemc_latitude,msldemc_longitude, ...
                        msldemc_imFOVmask_img,S_im,L_im,cmmdl_geo,K_L,K_S,1); toc;
                case 'L2PBK_LL0_M3'
                    tic; [ msldemc_imUFOVmask_ctr] =  ...
                        iaumars_get_msldemtUFOVmask_ctr_L2PBK_LL0_M3_mex(...
                        MSLDEMdata.imgpath,MSLDEMdata.hdr, MSLDEMdata.OFFSET, ...
                        msldemc_imFOVmask.chdr,...
                        msldemc_latitude,msldemc_longitude, ...
                        msldemc_imFOVmask_img,S_im,L_im,cmmdl_geo,K_L,K_S,0); toc;
                case 'L2PBK_LL0DYU_M3'
                    tic; [ msldemc_imUFOVmask_ctr] =  ...
                        iaumars_get_msldemtUFOVmask_ctr_L2PBK_LL0_M3_mex(...
                        MSLDEMdata.imgpath,MSLDEMdata.hdr, MSLDEMdata.OFFSET, ...
                        msldemc_imFOVmask.chdr,...
                        msldemc_latitude,msldemc_longitude, ...
                        msldemc_imFOVmask_img,S_im,L_im,cmmdl_geo,K_L,K_S,1); toc;
                case 'L2PBK_LLADYU_M3'
                    tic; [ msldemc_imUFOVmask_ctr] =  ...
                        iaumars_get_msldemtUFOVmask_ctr_L2PBK_LLA_M3_mex(...
                        MSLDEMdata.imgpath,MSLDEMdata.hdr, MSLDEMdata.OFFSET, ...
                        msldemc_imFOVmask.chdr,...
                        msldemc_latitude,msldemc_longitude, ...
                        msldemc_imFOVmask_img,S_im,L_im,cmmdl_geo,K_L,K_S); toc;
                case 'L2SMP_DAR_M3'
                    tic; [ msldemc_imUFOVmask_ctr] =  ...
                        iaumars_get_msldemtUFOVmask_ctr_L2SMP_DAR_M3_mex(...
                        MSLDEMdata.imgpath,MSLDEMdata.hdr, MSLDEMdata.OFFSET, ...
                        msldemc_imFOVmask.chdr,...
                        msldemc_latitude,msldemc_longitude, ...
                        msldemc_imFOVmask_img,S_im,L_im,cmmdl_geo,0); toc;
                case 'L2SMP_DARDYM_M3'
                    tic; [ msldemc_imUFOVmask_ctr] =  ...
                        iaumars_get_msldemtUFOVmask_ctr_L2SMP_DAR_M3_mex(...
                        MSLDEMdata.imgpath,MSLDEMdata.hdr, MSLDEMdata.OFFSET, ...
                        msldemc_imFOVmask.chdr,...
                        msldemc_latitude,msldemc_longitude, ...
                        msldemc_imFOVmask_img,S_im,L_im,cmmdl_geo,1); toc;
                case 'L2SMP_LL0_M3'
                    tic; [ msldemc_imUFOVmask_ctr] =  ...
                        iaumars_get_msldemtUFOVmask_ctr_L2SMP_LL0_M3_mex(...
                        MSLDEMdata.imgpath,MSLDEMdata.hdr, MSLDEMdata.OFFSET, ...
                        msldemc_imFOVmask.chdr,...
                        msldemc_latitude,msldemc_longitude, ...
                        msldemc_imFOVmask_img,S_im,L_im,cmmdl_geo,0); toc;
                case 'L2SMP_LL0DYU_M3'
                    tic; [ msldemc_imUFOVmask_ctr] =  ...
                        iaumars_get_msldemtUFOVmask_ctr_L2SMP_LL0_M3_mex(...
                        MSLDEMdata.imgpath,MSLDEMdata.hdr, MSLDEMdata.OFFSET, ...
                        msldemc_imFOVmask.chdr,...
                        msldemc_latitude,msldemc_longitude, ...
                        msldemc_imFOVmask_img,S_im,L_im,cmmdl_geo,1); toc;
                case 'L2SMP_LLP_M3_SCIMXY'
                    tic; [ msldemc_imUFOVmask_ctr] =  ...
                        iaumars_get_msldemtUFOVmask_ctr_L2SMP_LLP_M3_SCIMXY_mex(...
                        MSLDEMdata.imgpath,MSLDEMdata.hdr, MSLDEMdata.OFFSET, ...
                        msldemc_imFOVmask.chdr,...
                        msldemc_latitude,msldemc_longitude, ...
                        msldemc_imFOVmask_img,S_im,L_im,cmmdl_geo,0); toc;
                case 'L2SMP_LLPDYU_M3_SCIMXY'
                    tic; [ msldemc_imUFOVmask_ctr] =  ...
                        iaumars_get_msldemtUFOVmask_ctr_L2SMP_LLP_M3_SCIMXY_mex(...
                        MSLDEMdata.imgpath,MSLDEMdata.hdr, MSLDEMdata.OFFSET, ...
                        msldemc_imFOVmask.chdr,...
                        msldemc_latitude,msldemc_longitude, ...
                        msldemc_imFOVmask_img,S_im,L_im,cmmdl_geo,1); toc;
                case 'L2SMP_LLPDYU_M3_SCIMXY_STAA'
                    tic; [ msldemc_imUFOVmask_ctr] =  ...
                        iaumars_get_msldemtUFOVmask_ctr_L2SMP_LLP_M3_SCIMXY_STAA_mex(...
                        MSLDEMdata.imgpath,MSLDEMdata.hdr, MSLDEMdata.OFFSET, ...
                        msldemc_imFOVmask.chdr,...
                        msldemc_latitude,msldemc_longitude, ...
                        msldemc_imFOVmask_img,S_im,L_im,cmmdl_geo,1); toc;
                case 'L2SMP_LLP_M2_SCIMXY'
                    tic; [ msldemc_imUFOVmask_ctr] =  ...
                        iaumars_get_msldemtUFOVmask_ctr_L2SMP_LLP_M2_SCIMXY_mex(...
                        MSLDEMdata.imgpath,MSLDEMdata.hdr, MSLDEMdata.OFFSET, ...
                        msldemc_imFOVmask.chdr,...
                        msldemc_latitude,msldemc_longitude, ...
                        msldemc_imFOVmask_img,S_im,L_im,cmmdl_geo,0); toc;
                case 'L2SMP_LLPDYU_M2_SCIMXY'
                    tic; [ msldemc_imUFOVmask_ctr] =  ...
                        iaumars_get_msldemtUFOVmask_ctr_L2SMP_LLP_M2_SCIMXY_mex(...
                        MSLDEMdata.imgpath,MSLDEMdata.hdr, MSLDEMdata.OFFSET, ...
                        msldemc_imFOVmask.chdr,...
                        msldemc_latitude,msldemc_longitude, ...
                        msldemc_imFOVmask_img,S_im,L_im,cmmdl_geo,1); toc;
                
                otherwise
                    error('Undefined PROC_MODE %s',PROC_MODE);
            end
        otherwise
             error('COORDINATE_SYSTEM %s is not supported',coordsys);
    end
    
    if save_file
        data_type = envihdr_get_data_type_from_precision(class(msldemc_imUFOVmask_ctr));
        hdr_msldemc_imUFOVmask_ctr = mslgaleMosaicCrop_get_envihdr( ...
        MSLDEMdata,msldemc_imFOVmask.chdr,'BANDS',1,'DATA_TYPE',data_type);
        
        envi_save_raster(msldemc_imUFOVmask_ctr,hdr_msldemc_imUFOVmask_ctr, ...
            ufov_ctr_imgpath,ufov_ctr_hdrpath);
        
    end
else
    msldemc_imUFOVmask_ctr_obj = ENVIRasterSingleLayerMSLDEMCProj( ...
            basename_imUFOVmask_ctr,dirpath_cache);
    if ~strcmpi(msldemc_imUFOVmask_ctr_obj.hdr.msldem_basename,MSLDEMdata.basename)
        error('MSLDEM used for the cache is different from that of input.');
    end
    
    msldemc_imUFOVmask_ctr = msldemc_imUFOVmask_ctr_obj.readimg('precision','raw');
    
end    
msldemc_imUFOVmaskflg = or(msldemc_imUFOVmask_ctr,int8(msldemc_imFOVmask_pxlctrnn_img));

%% Safeguarding
%-------------------------------------------------------------------------%
% Evaluate the vertexes of each pixel
%-------------------------------------------------------------------------%
if contains(border_assess_opt,'d') % d means diagonal border
    basename_imUFOVmask_d = mastcam_get_basename_cache(basename_cache_com, ...
        'imUFOVmask_d',cache_vr);
    [ufov_d_imgpath,ufov_d_hdrpath] = mastcam_get_cachefilepath( ...
        basename_imUFOVmask_d,dirpath_cache);

    [flg_reproc] = doyouwanttoprocess({ufov_d_imgpath,ufov_d_hdrpath},force,load_cache_ifexist);
    
    if flg_reproc        
        % Evaluate the corner of the pixels
        %  o---------o
        %  |  (c,l)  |
        %  |    X    |  
        %  |         |
        %  o---------o (c+1/2,l+1/2)
        
        switch upper(coordsys)
            case {'NEE','NORTHEASTELEVATION','NORTHEASTNADIR','NENADIR'}
                msldemc_northing = MSLDEMdata.northing(lrange(1):lrange(2));
                msldemc_easting  = MSLDEMdata.easting(srange(1):srange(2));
                msldemc_xmc      = msldemc_northing - C_geo(1);
                msldemc_ymc      = msldemc_easting  - C_geo(2);
                tic; [ msldemt_imUFOVmask_d ] = get_msldemtUFOVmask_d_wmsldemc_L2_mex(...
                    MSLDEMdata.imgpath,MSLDEMdata.hdr,msldemc_imFOVmask.chdr,...
                    msldemc_xmc,msldemc_ymc,msldemc_imFOVmask_img,S_im,L_im,cmmdl_geo); toc;
            case {'IAU_MARS_SPHERE','IAU_MARS_ELLIPSOID'}
                msldemc_latitude  = deg2rad(MSLDEMdata.latitude(lrange(1):lrange(2)));
                msldemc_longitude = deg2rad(MSLDEMdata.longitude(srange(1):srange(2)));
                switch upper(PROC_MODE)
                    case 'L2PBK_LL0_M3'
                        tic; [ msldemt_imUFOVmask_d] =  ...
                            iaumars_get_msldemtUFOVmask_d_L2PBK_LL0_M3_mex(...
                            MSLDEMdata.imgpath,MSLDEMdata.hdr, MSLDEMdata.OFFSET, ...
                            msldemc_imFOVmask.chdr,...
                            msldemc_latitude,msldemc_longitude, ...
                            msldemc_imFOVmask_img,S_im,L_im,cmmdl_geo,K_L,K_S,0); toc;
                    case 'L2PBK_LL0DYU_M3'
                        tic; [ msldemt_imUFOVmask_d] =  ...
                            iaumars_get_msldemtUFOVmask_d_L2PBK_LL0_M3_mex(...
                            MSLDEMdata.imgpath,MSLDEMdata.hdr, MSLDEMdata.OFFSET, ...
                            msldemc_imFOVmask.chdr,...
                            msldemc_latitude,msldemc_longitude, ...
                            msldemc_imFOVmask_img,S_im,L_im,cmmdl_geo,K_L,K_S,1); toc;
                    otherwise
                        error('Undefined PROC_MODE %s',PROC_MODE);
                end
            otherwise
                error('COORDINATE_SYSTEM %s is not supported',coordsys);
        end
        
        tic; [ msldemc_imUFOVmask_d ] = ...
            get_msldemcUFOVmask_wVrtxPxlmsldemtUFOVmask(...
            msldemt_imUFOVmask_d,msldemc_imFOVmask_img); toc;
            
        if save_file
            data_type = envihdr_get_data_type_from_precision(class(msldemc_imUFOVmask_d));
            hdr_msldemc_imUFOVmask_d = mslgaleMosaicCrop_get_envihdr( ...
            MSLDEMdata,msldemc_imFOVmask.chdr,'BANDS',1,'DATA_TYPE',data_type);
        
            envi_save_raster(msldemc_imUFOVmask_d,hdr_msldemc_imUFOVmask_d, ...
                ufov_d_imgpath,ufov_d_hdrpath);
        end
    else
        msldemc_imUFOVmask_d_obj = ENVIRasterSingleLayerMSLDEMCProj( ...
        basename_imUFOVmask_d,dirpath_cache);
        if ~strcmpi(msldemc_imUFOVmask_d_obj.hdr.msldem_basename,MSLDEMdata.basename)
            error('MSLDEM used for the cache is different from that of input.');
        end

        msldemc_imUFOVmask_d = msldemc_imUFOVmask_d_obj.readimg('precision','raw');
    end
    msldemc_imUFOVmaskflg = int8(or(msldemc_imUFOVmask_d,msldemc_imUFOVmaskflg));
end

%%
%-------------------------------------------------------------------------%
% Evaluate the center of column edges
%-------------------------------------------------------------------------%
if contains(border_assess_opt,'c')
    error('future implementation');
%     basename_imUFOVmask_wcedge = sprintf('%s_imUFOVmask_wcedge_%s.mat',basename_cache,cache_vr);
%     fpath_imUFOVmask_wcedge = joinPath(dirpath_cache,basename_imUFOVmask_wcedge);
%     [flg_reproc] = doyouwanttoprocess(fpath_imUFOVmask_wcedge,force,load_cache_ifexist);
%     if flg_reproc
%         
%         % Evaluate the center of column edges
%         %   ---------
%         %  |  (c,l)  |
%         %  o    X    o
%         %  |         |
%         %   ---------
%         
%         tic; [dem_imxm,dem_imym] = get_imxycm_MSLDEM_mex(...
%             MSLDEMdata.imgpath,MSLDEMdata.hdr,MSTproj.msldemc_imFOVhdr,...
%             msldemc_northing+C_geo(1),msldemc_easting+C_geo(2),...
%             MSTproj.msldemc_imFOVmask,S_im,L_im,cmmdl_geo); toc;
% 
%         tic; [ msldemc_imUFOVmask_c ] = find_hidden_cedges_mastcamMSLDEM_v6_mex(...
%             msldemc_img,...0
%             msldemc_northing,...1
%             msldemc_easting,...2
%             MSTproj.msldemc_imFOVxy(:,:,1),...3
%             MSTproj.msldemc_imFOVxy(:,:,2),...4
%             MSTproj.msldemc_imFOVmask,...5
%             S_im,L_im,...6.7
%             dem_imxm,dem_imym); toc;
% 
%         clear dem_imxm dem_imym;
% 
%         % save msldemc_imUFOVmask_c.mat msldemc_imUFOVmask_c;
% 
%         tic; [ msldemc_imUFOVmask_wcedge ] = find_hidden_withcEdges_mastcamMSLDEM(...
%             msldemc_imUFOVmask_c,...0
%             MSTproj.msldemc_imFOVmask); toc; ...1
%             
%         if save_file
%             basename_msldem = MSLDEMdata.basename;
%             dirpath_msldem  = MSLDEMdata.dirpath;
%             if exist(fpath_imUFOVmask_wcedge,'file')
%                 delete(fpath_imUFOVmask_wcedge);
%             end
%             fprintf('Saving %s ...',fpath_imUFOVmask_wcedge);
%             save(fpath_imUFOVmask_wcedge,'msldemc_imUFOVmask_wcedge','msldemc_imFOVhdr',...
%                 'basename_msldem','dirpath_msldem');
%             fprintf('\nDone.\n');
%             
%         end
%     else
%         load(fpath_imUFOVmask_wcedge,'msldemc_imUFOVmask_wcedge','basename_msldem');
%         if ~strcmpi(basename_msldem,MSLDEMdata.basename)
%             error('MSLDEM used for the cache is different from that of input.');
%         end
%     end
%     msldemc_imUFOVmask = msldemc_imUFOVmask_wcedge + msldemc_imUFOVmask;
end
    
%%
%-------------------------------------------------------------------------%
% Evaluate the center of line edges
%-------------------------------------------------------------------------%
% Evaluate the center of line edges
        %   ----o----
        %  |  (c,l)  |
        %  |    X    |
        %  |         |
        %   ----o----
if contains(border_assess_opt,'l')
    error('future implementation');
%     basename_imUFOVmask_wledge = sprintf('%s_imUFOVmask_wledge_%s.mat',basename_cache,cache_vr);
%     fpath_imUFOVmask_wledge = joinPath(dirpath_cache,basename_imUFOVmask_wledge);
%     [flg_reproc] = doyouwanttoprocess(fpath_imUFOVmask_wledge,force,load_cache_ifexist);
%     if flg_reproc
%         
%         tic; [dem_imxm,dem_imym] = get_imxylm_MSLDEM_mex(...
%             MSLDEMdata.imgpath,MSLDEMdata.hdr,MSTproj.msldemc_imFOVhdr,...
%             msldemc_northing+C_geo(1),msldemc_easting+C_geo(2),...
%             MSTproj.msldemc_imFOVmask,S_im,L_im,cmmdl_geo); toc;
% 
% 
%         tic; [ msldemc_imUFOVmask_l ] = find_hidden_ledges_mastcamMSLDEM_v6_mex(...
%             msldemc_img,...0
%             msldemc_northing,...1
%             msldemc_easting,...2
%             MSTproj.msldemc_imFOVxy(:,:,1),...3
%             MSTproj.msldemc_imFOVxy(:,:,2),...4
%             MSTproj.msldemc_imFOVmask,...5
%             S_im,L_im,...6.7
%             dem_imxm,dem_imym); toc;
% 
%         clear dem_imxm dem_imym;
% 
%         % save msldemc_imUFOVmask_l.mat msldemc_imUFOVmask_l;
% 
%         tic; [ msldemc_imUFOVmask_wledge ] = find_hidden_withlEdges_mastcamMSLDEM(...
%             msldemc_imUFOVmask_l,...0
%             MSTproj.msldemc_imFOVmask); toc; ...1
%             
%         if save_file
%             basename_msldem = MSLDEMdata.basename;
%             dirpath_msldem  = MSLDEMdata.dirpath;
%             if exist(fpath_imUFOVmask_wledge,'file')
%                 delete(fpath_imUFOVmask_wledge);
%             end
%             fprintf('Saving %s ...',fpath_imUFOVmask_wledge);
%             save(fpath_imUFOVmask_wledge,'msldemc_imUFOVmask_wledge','msldemc_imFOVhdr',...
%                 'basename_msldem','dirpath_msldem');
%             fprintf('\nDone.\n');
%         end
%     else
%         load(fpath_imUFOVmask_wledge,'msldemc_imUFOVmask_wledge','basename_msldem');
%         if ~strcmpi(basename_msldem,MSLDEMdata.basename)
%             error('MSLDEM used for the cache is different from that of input.');
%         end
%     end
%     msldemc_imUFOVmask = msldemc_imUFOVmask_wledge + msldemc_imUFOVmask;
end

%% Summary
msldemc_imUFOVmask = msldemc_imFOVmask_img.*int8(msldemc_imUFOVmaskflg);
msldemc_imUFOVmask(msldemc_imUFOVmask==2) = 0;

end

