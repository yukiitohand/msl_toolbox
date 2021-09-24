function [msldemc_imUFOVxynn_obj] = mastcam_get_imUFOVxynn(mastcamdata_obj, ...
    MSLDEMdata,msldemc_imUFOVmask_obj,varargin)
% [msldemc_imUFOVxynn] = mastcam_get_imUFOVxynn(mastcamdata_obj, ...
%     MSLDEMdata,msldemc_imUFOVmask_obj,varargin)
%   Get msldemc_imUFOVmask filled with their rounded image coordinate.
%  INPUTS
%   mastcamdata_obj    : MASTCAMdata,MASTCAMgroup_wProcCode
%   MSLDEMdata         : MSLGaleDEMMosaic_v3
%   msldemc_imUFOVmask_obj : ENVIRasterSingleLayerMSLDEMCProj
%  OUTPUTS
%   msldemc_imUFOVxynn_obj: ENVIRasterMultBandMSLDEMCProj
%  OPTIONAL Parameters
%    'CACHE_DIRPATH': 
%       (default) msl_env_vars.dirpath_cache
%    'SAVE_FILE': 
%       (default) 1
%    'FORCE': boolean, whether to process forcefully
%       (default) 0
%    {'LOAD_CACHE_IFEXIST','LCIE'}
%       (default) 1
%    {'VARARGIN_PROCESS'}
%       cell array, varargin input for 
%             mastcam_calc_imUFOVxynn

global msl_env_vars
pdir_cache = msl_env_vars.dirpath_cache;
save_file = true;
force = false;
load_cache_ifexist = true;
cache_vr = ''; % {'v0','v1'}

varargin_proc = {};

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
            case {'VARARGIN_PROCESS'}
                varargin_proc = varargin{i+1};
                validateattributes(varargin_proc,{'cell'},{},mfilename,'VARARGIN_PROCESS');
            
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

if isempty(cache_vr)
   error('Please enter cache version with the option "CACHE_VER" or "CACHE_VERSION".'); 
end

% Input validation
validateattributes(mastcamdata_obj, ...
    {'MASTCAMdata','MASTCAMgroup_wProcCode'},{},mfilename,'mastcamdata_obj');
validateattributes(MSLDEMdata,{'MSLGaleDEMMosaic_v3'},{},mfilename,'MSLDEMdata');
validateattributes(msldemc_imUFOVmask_obj, ...
    {'ENVIRasterSingleLayerMSLDEMCProj','ENVIRasterMultBandMSLDEMCProj'}, ...
    {},mfilename,'msldemc_imUFOVmask_obj');

%%
[cache_dirname] = mastcam_get_cache_dirname(mastcamdata_obj);
dirpath_cache = joinPath(pdir_cache,cache_dirname);
if ~exist(dirpath_cache,'dir')
    mkdir(dirpath_cache);
end
[basename_cache_com] = mastcam_create_basename_cache(mastcamdata_obj);

basename_imUFOVxynn = mastcam_get_basename_cache(basename_cache_com, ...
    'imUFOVxynn',cache_vr);
[imUFOVxynn_imgpath,imUFOVxynn_hdrpath] = mastcam_get_cachefilepath( ...
    basename_imUFOVxynn,dirpath_cache);

% Evaluate to perform the 
[flg_reproc] = doyouwanttoprocess({imUFOVxynn_imgpath,imUFOVxynn_hdrpath}, ...
    force,load_cache_ifexist);

if flg_reproc
    % Mask invisible pixels
    [msldemc_imUFOVxynn_img] = mastcam_calc_imUFOVxynn(mastcamdata_obj, ...
        MSLDEMdata, msldemc_imUFOVmask_obj,varargin_proc{:});
    
    [data_type] = envihdr_get_data_type_from_precision(class(msldemc_imUFOVxynn_img));
    hdr_imUFOVxynn = mslgaleMosaicCrop_get_envihdr(MSLDEMdata, ...
        msldemc_imUFOVmask_obj.chdr,'BANDS',2,'DATA_TYPE',data_type);
    
    if save_file
        envi_save_raster(msldemc_imUFOVxynn_img,hdr_imUFOVxynn, ...
            imUFOVxynn_imgpath,imUFOVxynn_hdrpath);
        msldemc_imUFOVxynn_obj = ENVIRasterMultBandMSLDEMCProj( ...
            basename_imUFOVxynn,dirpath_cache);
    else
        msldemc_imUFOVxynn_obj = ENVIRasterMultBandMSLDEMCProj('','');
        msldemc_imUFOVxynn_obj.img = msldemc_imUFOVxynn_img;
        msldemc_imUFOVxynn_obj.hdr = hdr_imUFOVxynn;
        msldemc_imUFOVxynn_obj.chdr = msldemc_imUFOVhdr;
        msldemc_imUFOVxynn_obj.proj_info = envi_get_proj_info_SphereEquirect( ...
            hdr_imUFOVxynn);
    end
    
else
    msldemc_imUFOVxynn_obj = ENVIRasterMultBandMSLDEMCProj( ...
        basename_imUFOVxynn,dirpath_cache);
    
    if ~strcmpi(msldemc_imUFOVxynn_obj.hdr.msldem_basename,MSLDEMdata.basename)
        error('MSLDEM used for the cache is different from that of input.');
    end
end

end