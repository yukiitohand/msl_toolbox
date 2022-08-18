function [msldemc_imUFOVmask] = mastcam_getUFOVwMSLDEM_traditional( ...
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
%      msldemc_imFOVmask
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
%   "PROC_MODE": char array, processing mode. {'XDraw','R3','R2'}
%     (default) 'XDRAW'

global msl_env_vars
pdir_cache = msl_env_vars.dirpath_cache;
save_file = true;
force = false;
load_cache_ifexist = true;
cache_vr = ''; % {'v0','v1'}
PROC_MODE = 'XDraw';

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
            case 'PROC_MODE'
                PROC_MODE = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

if isempty(cache_vr)
   error('Please enter cache version with the option "CACHE_VER" or "CACHE_VERSION".'); 
end

%-------------------------------------------------------------------------%
% Get cam_C_geo
%-------------------------------------------------------------------------%
cmmdl = mastcamdata_obj.CAM_MDL;
rover_nav_coord = mastcamdata_obj.ROVER_NAV;
cmmdl_geo = transform_CAHVOR_MODEL_wROVER_NAV(cmmdl,rover_nav_coord);

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

%-------------------------------------------------------------------------%
% convert camera center to the MSLDEMC image coordinate system
%-------------------------------------------------------------------------%
cmmdl_geo.C(1) = msldemc_imFOVmask.northing2y(cmmdl_geo.C(1));
cmmdl_geo.C(2) = msldemc_imFOVmask.easting2x(cmmdl_geo.C(2));
cmmdl_geo.C(3) = -cmmdl_geo.C(3);

%%
%-------------------------------------------------------------------------%
% Evaluate the center of each pixel
%-------------------------------------------------------------------------%
basename_imUFOVmask_ctr = mastcam_get_basename_cache(basename_cache_com, ...
    sprintf('imUFOVmask_%s_ctr',PROC_MODE),cache_vr);
[ufov_ctr_imgpath,ufov_ctr_hdrpath] = mastcam_get_cachefilepath( ...
    basename_imUFOVmask_ctr,dirpath_cache);

[flg_reproc] = doyouwanttoprocess({ufov_ctr_imgpath,ufov_ctr_hdrpath},force,load_cache_ifexist);
if flg_reproc
    
            
    switch upper(PROC_MODE)
        case 'XDRAW'
            tic; [ msldemc_imUFOVmask_ctr ] = msldemc_get_viewshed_XDraw_mex(...
                MSLDEMdata.imgpath,MSLDEMdata.hdr,msldemc_imFOVmask.chdr,...
                msldemc_imFOVmask_img,cmmdl_geo); toc;
        case 'R3'
            srange = msldemc_imFOVmask.get_xrange_base();
            lrange = msldemc_imFOVmask.get_yrange_base();
            msldemc_img = MSLDEMdata.get_subimage_wPixelRange(srange,lrange,'precision','double');
            tic; [ msldemc_imUFOVmask_ctr ] = msldemc_get_viewshed_R3_mex(...
                msldemc_img,msldemc_imFOVmask_img,cmmdl_geo); toc;

        otherwise
            error('Undefined PROC_MODE %s',PROC_MODE);
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




%% Summary
msldemc_imUFOVmask = msldemc_imFOVmask_img.*int8(msldemc_imUFOVmaskflg);
msldemc_imUFOVmask(msldemc_imUFOVmask==2) = 0;

end

