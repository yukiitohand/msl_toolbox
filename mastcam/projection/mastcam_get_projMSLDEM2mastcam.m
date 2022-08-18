function [msldemc_imFOVmask,msldemc_imFOVxy,msldemc_imFOVhdr]...
    = mastcam_get_projMSLDEM2mastcam(MSLDEMdata,mstdata_obj,varargin)

global msl_env_vars
cachedirpath = msl_env_vars.dirpath_cache;
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
                cachedirpath = varargin{i+1};
            case 'SAVE_FILE'
                save_file = varargin{i+1};
            case 'FORCE'
                force = varargin{i+1};
            case {'LOAD_CACHE_IFEXIST','LCIE'}
                load_cache_ifexist = varargin{i+1};
            case {'CACHE_VER','CACHE_VERSION'}
                cache_vr = varargin{i+1};
                
            % ## PROCESSING OPTIONS #--------------------------------------
            case 'VARARGIN_PROCESS'
                varargin_proc = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

if isempty(cache_vr)
   error('Please enter cache version with the option "CACHE_VER" or "CACHE_VERSION".'); 
end

[basename_cache_com] = mastcam_create_basename_cache(mstdata_obj);
cachepath = joinPath(cachedirpath,...
    [basename_cache_com sprintf('_projMSLDEM2mastcam_%s.mat',cache_vr)]);

[flg_reproc] = doyouwanttoprocess(cachepath,force,load_cache_ifexist);

if flg_reproc
    [msldemc_imFOVmask,msldemc_imFOVhdr,L_im,S_im,cmmdl_geo,option]...
        = mastcam_get_projMSLDEM2mastcam_v3_imFOVmask(MSLDEMdata,...
        mstdata_obj,varargin_proc{:});
    
    if save_file
        basename_msldem = MSLDEMdata.basename;
        dirpath_msldem  = MSLDEMdata.dirpath;
        if exist(cachepath,'file')
            delete(cachepath);
        end
        fprintf('Saving %s ...',cachepath);
        save(cachepath,'msldemc_imFOVmask','msldemc_imFOVhdr',...
            'L_im','S_im','cmmdl_geo','option',...
            'basename_msldem','dirpath_msldem',... 
            '-nocompression');
        fprintf('\nDone.\n');
    end
    
else
    load(cachepath,'msldemc_imFOVmask','msldemc_imFOVhdr','L_im','S_im',...
        'cmmdl_geo','basename_msldem','dirpath_msldem');
    
    if ~strcmpi(basename_msldem,MSLDEMdata.basename)
        error('MSLDEM used for the cache is different from that of input.');
    end
end

% imxy is no longer calculated here because it makes the later processing
% slower.
msldemc_imFOVxy = [];

% % does not save imxy since it can be too big, loading takes much longer.
% switch mstdata_obj.Linearization
%     case 1
%         [msldemc_imFOVxy] = mastcam_get_projMSLDEM2mastcam_v3_imxy(...
%             MSLDEMdata,msldemc_imFOVmask,msldemc_imFOVhdr,L_im,S_im,cmmdl_geo);
%     case 0
%         msldemc_imFOVxy = [];
%     otherwise
%         error('Not defined for Linearization=%d',mstdata_obj.Linearization);
% end

end



