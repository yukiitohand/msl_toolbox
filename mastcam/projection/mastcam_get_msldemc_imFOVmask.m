function [msldemc_imFOVmask_obj,suppl_info] = mastcam_get_msldemc_imFOVmask( ...
    MSLDEMdata,mstdata_obj,varargin)
% [msldemc_imFOVmask_obj,suppl_info] = mastcam_get_msldemc_imFOVmask( ...
%     MSLDEMdata,mstdata_obj,varargin)
% Wrapper function for getting msldemc_imFOVmask.
%  Actual computation is performed by 
%    mastcam_get_projMSLDEM2mastcam_v3_imFOVmask
% Wrapper function reads the data from the saved file if exists by default.
%
%  INPUTS
%    MSLDEMdata: MSLGaleDEMMosaic_v3 obj, DEM data, which which FOVmask is
%       calculated.
%    mstdata_obj: MASTCAMdata or MASTCAMgroup_eye object for which
%       CAHV/CAHVOR model is same.
%  OUTPUTS
%    msldemc_imFOVmask_obj: ENVIRasterSingleLayerMSLDEMCProj obj
%    suppl_info: struct, supplemental information
%       'L_im','S_im','cmmdl_geo','option'
%  REQUIRED Parameters
%    {'CACHE_VER','CACHE_VERSION'}
%        char, cache version, like 'v0','v1'
%  OPTIONAL Parameters
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
%            mastcam_get_projMSLDEM2mastcam_v3_imFOVmask
%       (default) {}
%
% (c) 2021 Yuki Itoh.


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
            case 'VARARGIN_PROCESS'
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
validateattributes(mstdata_obj, ...
    {'MASTCAMdata','MASTCAMgroup_wProcCode'},{},mfilename,'mstdata_obj');
validateattributes(MSLDEMdata,{'MSLGaleDEMMosaic_v3'},{},mfilename,'MSLDEMdata');

%%
[cache_dirname] = mastcam_get_cache_dirname(mstdata_obj);
dirpath_cache = joinPath(pdir_cache,cache_dirname);
if ~exist(dirpath_cache,'dir')
    mkdir(dirpath_cache);
end
[basename_cache_com] = mastcam_create_basename_cache(mstdata_obj);
basename_imFOVmask = mastcam_get_basename_cache(basename_cache_com,'imFOVmask',cache_vr);
[imgpath,hdrpath] = mastcam_get_cachefilepath(basename_imFOVmask,dirpath_cache);
supplpath = joinPath(dirpath_cache,[basename_imFOVmask '_suppl.mat']);

[flg_reproc] = doyouwanttoprocess({imgpath,hdrpath,supplpath},force,load_cache_ifexist);

if flg_reproc
    [msldemc_imFOVmask,msldemc_imFOVhdr,L_im,S_im,cmmdl_geo,option]...
        = mastcam_get_projMSLDEM2mastcam_v3_imFOVmask(MSLDEMdata,...
        mstdata_obj,varargin_proc{:});
    
    [data_type] = envihdr_get_data_type_from_precision('int8');
    hdr_c = mslgaleMosaicCrop_get_envihdr(MSLDEMdata,msldemc_imFOVhdr, ...
        'BANDS',1,'DATA_TYPE',data_type);
    
    if nargout>1
        suppl_info = struct('L_im',L_im,'S_im',S_im, ...
            'cmmdl_geo',cmmdl_geo,'option',option);
    end
    
    if save_file
        envi_save_raster(msldemc_imFOVmask,hdr_c,imgpath,hdrpath);
        
        if exist(supplpath,'file'),delete(supplpath); end
        fprintf('Saving %s ...',supplpath);
        save(supplpath,'L_im','S_im','cmmdl_geo','option');
        fprintf('\nDone.\n');
        
        msldemc_imFOVmask_obj = ENVIRasterSingleLayerMSLDEMCProj( ...
            basename_imFOVmask,dirpath_cache);
        
    else
        msldemc_imFOVmask_obj = ENVIRasterSingleLayerMSLDEMCProj('','');
        msldemc_imFOVmask_obj.img = msldemc_imFOVmask;
        msldemc_imFOVmask_obj.hdr = hdr_c;
        msldemc_imFOVmask_obj.chdr = msldemc_imFOVhdr;
        msldemc_imFOVmask_obj.proj_info = envi_get_proj_info_SphereEquirect( ...
            msldemc_imFOVmask_obj.hdr);
    end
    
else
    msldemc_imFOVmask_obj = ENVIRasterSingleLayerMSLDEMCProj(basename_imFOVmask,dirpath_cache);
    
    if ~strcmpi(msldemc_imFOVmask_obj.hdr.msldem_basename,MSLDEMdata.basename)
        error('MSLDEM used for the cache is different from that of input.');
    end
    
    if nargout>1
        suppl_info = load(supplpath);
    end
    
end

end



