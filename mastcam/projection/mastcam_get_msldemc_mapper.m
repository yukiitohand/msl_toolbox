function [mastcam_msldemc_mapper_obj] = mastcam_get_msldemc_mapper( ...
    mstdata_obj,MSLDEMdata,msldemc_imUFOVxynn_obj,mastcam_nn_msldem,varargin)
% [mastcam_msldemc_mapper_obj] = mastcam_get_msldemc_mapper( ...
%     mstdata_obj,MSLDEMdata,msldemc_imUFOVxynn_obj,mastcam_nn_msldem,varargin)
%  Create a mastcam <-> msldem two-way map projection.
% INPUTS
%   mstdata_obj: either of 'MASTCAMdata' or 'MASTCAMgroup_wProcCode'
%   MSLDEMdata : MSLDEMGaleMosaic_v3
%   msldemc_imUFOVxynn_obj : ENVIRasterMultBandMSLDEMCProj
%   mastcam_nn_msldem : ENVIRasterMultBand
% OUTPUTS
%   mastcam_msldemc_mapper_obj: MASTCAM_MSLDEMC_Mapper class
% OPTIONAL Parameters
%    'CACHE_DIRPATH': 
%       (default) msl_env_vars.dirpath_cache
%    'SAVE_FILE': 
%       (default) 1
%    'FORCE': boolean, whether to process forcefully
%       (default) 0
%    {'LOAD_CACHE_IFEXIST','LCIE'}
%       (default) 1
%
%

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
            case {'VARARGIN_PROCESS'}
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

if isempty(cache_vr)
   error('Please enter cache version with the option "CACHE_VER" or "CACHE_VERSION".'); 
end

%
validateattributes(mstdata_obj, ...
    {'MASTCAMdata','MASTCAMgroup_wProcCode'},{},mfilename,'mstdata_obj');
validateattributes(MSLDEMdata,{'MSLGaleDEMMosaic_v3'},{},mfilename,'MSLDEMdata');
validateattributes(msldemc_imUFOVxynn_obj, ...
    {'ENVIRasterMultBandMSLDEMCProj'},{},mfilename,'msldemc_imUFOVxynn_obj');
validateattributes(mastcam_nn_msldem, ...
    {'ENVIRasterMultBand'},{},mfilename,'mastcam_nn_msldem');

%%
[cache_dirname] = mastcam_get_cache_dirname(mstdata_obj);
dirpath_cache = joinPath(pdir_cache,cache_dirname);
if ~exist(dirpath_cache,'dir')
    mkdir(dirpath_cache);
end
[basename_cache_com] = mastcam_create_basename_cache(mstdata_obj);
basename_msldemc_mapper = mastcam_get_basename_cache(basename_cache_com, ...
    'msldemc_mapper',cache_vr);
mastcam_msldemc_mapper_filepath = joinPath(dirpath_cache, ...
    [basename_msldemc_mapper '.mat']);


% Evaluate to perform the 
[flg_reproc] = doyouwanttoprocess({mastcam_msldemc_mapper_filepath}, ...
    force,load_cache_ifexist);

if flg_reproc
    % Mask invisible pixels
    [mapper_mastcam2msldemc,mapper_msldemc2mastcam_mat, ...
        mapcell_msldemc2mastcam] = mastcam_create_msldemc_mapper( ...
        msldemc_imUFOVxynn_obj,mastcam_nn_msldem);
    
    chdr = msldemc_imUFOVxynn_obj.chdr;
    proj_info = msldemc_get_proj_info_from_base(MSLDEMdata,msldemc_imUFOVxynn_obj);
    msldem_basename = msldemc_imUFOVxynn_obj.hdr.msldem_basename;
    msldem_dirpath  = msldemc_imUFOVxynn_obj.hdr.msldem_dirpath;
    if save_file
        save(mastcam_msldemc_mapper_filepath,'mapper_mastcam2msldemc', ...
            'mapper_msldemc2mastcam_mat','mapcell_msldemc2mastcam',...
            'chdr','proj_info','msldem_basename','msldem_dirpath');   
    end
else
    load(mastcam_msldemc_mapper_filepath,'mapper_mastcam2msldemc', ...
            'mapper_msldemc2mastcam_mat','mapcell_msldemc2mastcam', ...
            'chdr','proj_info','msldem_basename','msldem_dirpath');
    if ~strcmpi(msldem_basename,MSLDEMdata.basename)
        error('MSLDEM used for the cache is different from that of input.');
    end
end

mastcam_msldemc_mapper_obj = MASTCAM_MSLDEMC_Mapper();
mastcam_msldemc_mapper_obj.mastcam2msldemc = mapper_mastcam2msldemc;
mastcam_msldemc_mapper_obj.msldemc2mastcam_mat = mapper_msldemc2mastcam_mat;
mastcam_msldemc_mapper_obj.mapcell_msldemc2mastcam = mapcell_msldemc2mastcam;
mastcam_msldemc_mapper_obj.msldem_basename = msldem_basename;
mastcam_msldemc_mapper_obj.msldem_dirpath  = msldem_dirpath;
mastcam_msldemc_mapper_obj.msldemc_chdr = chdr;
mastcam_msldemc_mapper_obj.msldemc_proj_info = proj_info;


%%


end