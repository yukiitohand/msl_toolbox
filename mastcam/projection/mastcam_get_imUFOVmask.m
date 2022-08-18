function [msldemc_imUFOVmask_obj] = mastcam_get_imUFOVmask( ...
    mastcamdata_obj,MSLDEMdata,msldemc_imFOVmask_obj,msldemc_imFOVmask_pxlctrnn_obj,varargin)
% [msldemc_imUFOVmask_obj] = mastcam_get_imUFOVmask( ...
%     mastcamdata_obj,MSLDEMdata,msldemc_imFOVmask,msldemc_imFOVmask_pxlctrnn_obj,varargin)
%  Mask invisible points within imFOVmask and crop the marginal pixels.
% INPUTS
%   mastcamdata_obj: either of 'MASTCAMdata' or 'MASTCAMgroup_wProcCode'
%   MSLDEMdata : MSLDEMGaleMosaic_v3
%   msldemc_imFOVmask_obj : ENVIRasterSingleLayerMSLDEMCProj
%   msldemc_imFOVmask_pxlctrnn_obj: ENVIRasterSingleLayerMSLDEMCProj
% OUTPUTS
%   msldemc_imUFOVmask_obj : ENVIRasterSingleLayerMSLDEMCProj
% 
% OPTIONAL Parameters
%    'CACHE_DIRPATH': 
%       (default) msl_env_vars.dirpath_cache
%    'SAVE_FILE': 
%       (default) 1
%    'FORCE': boolean, whether to process forcefully
%       (default) 0
%    {'LOAD_CACHE_IFEXIST','LCIE'}
%       (default) 1
%    'VARARGIN_PROCESS'
%       cell array, Optional Parameters for 
%            proj_mastcam2MSLDEM_v5_mexw_v2
%       (default) {}
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
save_file = true;
force = false;
load_cache_ifexist = true;
cache_vr = ''; % {'v0','v1'}

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
                validateattributes(save_file,{'numeric','logical'}, ...
                    {'binary'},mfilename,'SAVE_FILE');
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
            case {'BORDER_ASSESS_OPT','COORDINATE_SYSTEM','PROC_MODE','K_L','K_S'}
                
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

if isempty(cache_vr)
   error('Please enter cache version with the option "CACHE_VER" or "CACHE_VERSION".'); 
end

% INPUT validation
validateattributes(mastcamdata_obj, ...
    {'MASTCAMdata','MASTCAMgroup_wProcCode'},{},mfilename,'mastcamdata_obj');
validateattributes(MSLDEMdata,{'MSLGaleDEMMosaic_v3'},{},mfilename,'MSLDEMdata');
validateattributes(msldemc_imFOVmask_obj, ...
    {'ENVIRasterSingleLayerMSLDEMCProj','ENVIRasterMultBandMSLDEMCProj'}, ...
    {},mfilename,'msldemc_imFOVmask');
validateattributes(msldemc_imFOVmask_pxlctrnn_obj, ...
    {'ENVIRasterSingleLayerMSLDEMCProj','ENVIRasterMultBandMSLDEMCProj'}, ...
    {},mfilename,'msldemc_imFOVmask_pxlctrnn_obj');


%%
[cache_dirname] = mastcam_get_cache_dirname(mastcamdata_obj);
dirpath_cache = joinPath(pdir_cache,cache_dirname);
if ~exist(dirpath_cache,'dir')
    mkdir(dirpath_cache);
end
[basename_cache_com] = mastcam_create_basename_cache(mastcamdata_obj);

basename_imUFOVmask = mastcam_get_basename_cache(basename_cache_com,'imUFOVmask',cache_vr);
[imUFOVmask_imgpath,imUFOVmask_hdrpath] = mastcam_get_cachefilepath( ...
    basename_imUFOVmask,dirpath_cache);

% Evaluate to perform the 
[flg_reproc] = doyouwanttoprocess({imUFOVmask_imgpath,imUFOVmask_hdrpath}, ...
    force,load_cache_ifexist);

if flg_reproc
    % Mask invisible pixels
    [msldemc_imUFOVmask] = mastcam_getUFOVwMSLDEM_mexw_v2( ...
        mastcamdata_obj,MSLDEMdata,msldemc_imFOVmask_obj, ...
        msldemc_imFOVmask_pxlctrnn_obj,varargin{:});
    
    % crop the mask
    [msldemc_imUFOVmask_crop,msldemc_imUFOVhdr] ...
        = mastcam_crop_msldemc_imUFOVmask_v2( ...
        msldemc_imFOVmask_obj,msldemc_imUFOVmask);
    
    [data_type] = envihdr_get_data_type_from_precision('int8');
    hdr_imUFOVmask = mslgaleMosaicCrop_get_envihdr(MSLDEMdata, ...
        msldemc_imUFOVhdr, 'BANDS',1,'DATA_TYPE',data_type);
    
    if save_file
        envi_save_raster(msldemc_imUFOVmask_crop,hdr_imUFOVmask,imUFOVmask_imgpath,imUFOVmask_hdrpath);
        msldemc_imUFOVmask_obj = ENVIRasterSingleLayerMSLDEMCProj( ...
            basename_imUFOVmask,dirpath_cache);
    else
        msldemc_imUFOVmask_obj = ENVIRasterSingleLayerMSLDEMCProj('','');
        msldemc_imUFOVmask_obj.img = msldemc_imUFOVmask_crop;
        msldemc_imUFOVmask_obj.hdr = hdr_imUFOVmask;
        msldemc_imUFOVmask_obj.chdr = msldemc_imUFOVhdr;
        msldemc_imUFOVmask_obj.proj_info = envi_get_proj_info_SphereEquirect( ...
            hdr_imUFOVmask);
    end
    
else
    msldemc_imUFOVmask_obj = ENVIRasterSingleLayerMSLDEMCProj( ...
        basename_imUFOVmask,dirpath_cache);
    
    if ~strcmpi(msldemc_imUFOVmask_obj.hdr.msldem_basename,MSLDEMdata.basename)
        error('MSLDEM used for the cache is different from that of input.');
    end
end

end