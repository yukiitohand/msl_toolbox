function [msldemc_imFOVmaskr13_obj] = mastcam_get_msldemc_imFOVmask_13refined( ...
    mstdata_obj,msldemc_imFOVmask_obj,mastcam_ref_msldem_obj,varargin)
% [msldemc_imFOVmaskr13_obj] = mastcam_get_msldemc_imFOVmask_13refined( ...
%     mstdata_obj,msldemc_imFOVmask_obj,mastcam_ref_msldemc_obj,varargin)
% Evaluate the pixels in the msldemc_imFOVmask that are marked as 3.
%  INPUTS
%    mstdata_obj: MASTCAMdata or MASTCAMgroup_wProcCode class object
%    msldemc_imFOVmask_obj: ENVIRasterSingleLayerMSLDEMCProj class object
%    mastcam_ref_msldem_obj: ENVIRasterMultBand class object
%  OUTPUTS
%    msldemc_imFOVmaskr13_obj: ENVIRasterSingleLayerMSLDEMCProj class object
%
%  OPTIONAL Parameters
%    'CACHE_DIRPATH': 
%       (default) msl_env_vars.dirpath_cache
%    'SAVE_FILE': 
%       (default) 1
%    'FORCE': boolean, whether to process forcefully
%       (default) 0
%    {'LOAD_CACHE_IFEXIST','LCIE'}
%       (default) 1

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
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

if isempty(cache_vr)
   error('Please enter cache version with the option "CACHE_VER" or "CACHE_VERSION".'); 
end

% INPUT validation
validateattributes(mstdata_obj, ...
    {'MASTCAMdata','MASTCAMgroup_wProcCode'},{},mfilename,'mstdata_obj');
validateattributes(msldemc_imFOVmask_obj, ...
    {'ENVIRasterSingleLayerMSLDEMCProj','ENVIRasterMultBandMSLDEMCProj'}, ...
    {},mfilename,'msldemc_imFOVmask');
validateattributes(mastcam_ref_msldem_obj, ...
    {'ENVIRasterMultBand'},{},mfilename,'mastcam_ref_msldem_obj');

%%
[cache_dirname] = mastcam_get_cache_dirname(mstdata_obj);
dirpath_cache   = joinPath(pdir_cache,cache_dirname);
if ~exist(dirpath_cache,'dir')
    mkdir(dirpath_cache);
end
[basename_cache_com]  = mastcam_create_basename_cache(mstdata_obj);
basename_imFOVmaskr13 = mastcam_get_basename_cache(basename_cache_com,'imFOVmaskr13',cache_vr);
[imgpath,hdrpath] = mastcam_get_cachefilepath(basename_imFOVmaskr13,dirpath_cache);

[flg_reproc] = doyouwanttoprocess({imgpath,hdrpath},force,load_cache_ifexist);
%%
if flg_reproc
    msldemc_imFOVmaskr13_img = mastcam_msldemc_imFOVmask_eval_13_v2( ...
        msldemc_imFOVmask_obj,mastcam_ref_msldem_obj);
    
    hdr_r13 = msldemc_imFOVmask_obj.hdr;
    
    if save_file
        envi_save_raster(msldemc_imFOVmaskr13_img,hdr_r13,imgpath,hdrpath);
        msldemc_imFOVmaskr13_obj = ENVIRasterSingleLayerMSLDEMCProj( ...
            basename_imFOVmaskr13,dirpath_cache);
        
    else
        msldemc_imFOVmaskr13_obj = ENVIRasterSingleLayerMSLDEMCProj('','');
        msldemc_imFOVmaskr13_obj.img  = msldemc_imFOVmaskr13_img;
        msldemc_imFOVmaskr13_obj.hdr  = hdr_r13;
        msldemc_imFOVmaskr13_obj.chdr = msldemc_imFOVmask_obj.chdr;
        msldemc_imFOVmaskr13_obj.proj_info = envi_get_proj_info_SphereEquirect( ...
            msldemc_imFOVmask_obj.hdr);
    end
    
else
    msldemc_imFOVmaskr13_obj = ENVIRasterSingleLayerMSLDEMCProj( ...
        basename_imFOVmaskr13,dirpath_cache);
    
    if ~strcmpi(msldemc_imFOVmaskr13_obj.hdr.msldem_basename, ...
            msldemc_imFOVmask_obj.hdr.msldem_basename)
        error('MSLDEM used for the cache is different from that of input.');
    end
    
end


end