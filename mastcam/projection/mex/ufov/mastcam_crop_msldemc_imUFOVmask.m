function [msldemc_imUFOVmask,msldemc_imUFOVhdr,msldemc_imUFOVxynn,...
    mastcam_msldemc_nn_UFOVmask] = mastcam_crop_msldemc_imUFOVmask(MSTproj,varargin)

global msl_env_vars
dirpath_cache = msl_env_vars.dirpath_cache;
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
            case {'CACHE_DIRPATH','DIRPATH_CACHE'}
                dirpath_cache = varargin{i+1};
            case 'SAVE_FILE'
                save_file = varargin{i+1};
            case 'FORCE'
                force = varargin{i+1};
            case {'LOAD_CACHE_IFEXIST','LCIE'}
                load_cache_ifexist = varargin{i+1};
            case {'CACHE_VER','CACHE_VERSION'}
                cache_vr = varargin{i+1};
            
            % ## PROCESSING OPTIONS #--------------------------------------
            
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

if isempty(cache_vr)
   error('Please enter cache version with the option "CACHE_VER" or "CACHE_VERSION".'); 
end

mastcamdata_obj = MSTproj.MASTCAMdata;
MSLDEMdata = MSTproj.MSLDEMdata;
[basename_cache] = mastcam_create_basename_cache(mastcamdata_obj);

basename_cache_ufovc = sprintf('%s_imUFOV_%s.mat',basename_cache,cache_vr);
fpath_ufovc = joinPath(dirpath_cache,basename_cache_ufovc);
[flg_reproc] = doyouwanttoprocess(fpath_ufovc,force,load_cache_ifexist);

if flg_reproc
    if isempty(MSTproj.msldemc_imUFOVmask) || isempty(MSTproj.msldemc_imFOVhdr) ...
            || isempty(MSTproj.msldemc_imFOVxy) || isempty(MSTproj.mastcam_msldemc_nn)
        error('Processing required preceding computation is loaded');
    end
    
    % Detect the margin
    valid_lines   = find(any(MSTproj.msldemc_imUFOVmask',1));
    lrnge         = [valid_lines(1), valid_lines(end)];
    len_vl        = lrnge(2)-lrnge(1)+1;
    valid_samples = find(any(MSTproj.msldemc_imUFOVmask,1));
    srnge         = [valid_samples(1), valid_samples(end)];
    len_vs        = srnge(2)-srnge(1)+1;

    l1   = MSTproj.msldemc_imFOVhdr.line_offset+lrnge(1);
    lend = MSTproj.msldemc_imFOVhdr.line_offset+lrnge(2);
    s1   = MSTproj.msldemc_imFOVhdr.sample_offset+srnge(1);
    send = MSTproj.msldemc_imFOVhdr.sample_offset+srnge(2);
    msldemc_northing = MSLDEMdata.hdr.y(l1:lend);
    msldemc_easting  = MSLDEMdata.hdr.x(s1:send);

    % Crop the mask
    msldemc_imUFOVmask = MSTproj.msldemc_imUFOVmask(lrnge(1):lrnge(2),srnge(1):srnge(2));

    % Updated header information (including offset and the size of the mask
    % image
    msldemc_imUFOVhdr = [];
    msldemc_imUFOVhdr.lines   = len_vl;
    msldemc_imUFOVhdr.samples = len_vs;
    msldemc_imUFOVhdr.line_offset   = l1-1;
    msldemc_imUFOVhdr.sample_offset = s1-1;
    msldemc_imUFOVhdr.y = msldemc_northing;
    msldemc_imUFOVhdr.x = msldemc_easting;

    % Crop imxy, and fill NaN for the pixels that are not observed.
    msldemc_imUFOVxy = MSTproj.msldemc_imFOVxy(lrnge(1):lrnge(2),srnge(1):srnge(2),:);
    mm = (msldemc_imUFOVmask==0);
    for i=1:2
        msldemc_imUFOVtmp = msldemc_imUFOVxy(:,:,i);
        msldemc_imUFOVtmp(mm) = nan;
        msldemc_imUFOVxy(:,:,i) = msldemc_imUFOVtmp;
    end

    % Mapping from msldemc to the nearest pixel of the MASTCAM image
    msldemc_imUFOVxynn = round(msldemc_imUFOVxy+1);

    L_im = mastcamdata_obj.L_im; S_im = mastcamdata_obj.S_im;
    for i=1:2
        msldemc_imUFOVtmp = msldemc_imUFOVxynn(:,:,i);
        msldemc_imUFOVtmp(msldemc_imUFOVtmp<1) = 1;
        msldemc_imUFOVtmp(isnan(msldemc_imUFOVtmp)) = -1;
        if i==1
            msldemc_imUFOVtmp(msldemc_imUFOVtmp>S_im) = S_im;
        elseif i==2
            msldemc_imUFOVtmp(msldemc_imUFOVtmp>L_im) = L_im;
        end
        msldemc_imUFOVxynn(:,:,i) = msldemc_imUFOVtmp;
    end
    msldemc_imUFOVxynn = int16(msldemc_imUFOVxynn);

    % pixel correspondence of mastcam_msldemc_nn is updated for the cropped
    % image.
    mastcam_msldemc_nn_UFOVmask = MSTproj.mastcam_msldemc_nn;
    mastcam_msldemc_nn_UFOVmask(:,:,1) = mastcam_msldemc_nn_UFOVmask(:,:,1)-int32(srnge(1)-1);
    mastcam_msldemc_nn_UFOVmask(:,:,2) = mastcam_msldemc_nn_UFOVmask(:,:,2)-int32(lrnge(1)-1);
    
    if save_file
        basename_msldem = MSLDEMdata.basename;
        dirpath_msldem  = MSLDEMdata.dirpath;
        if exist(fpath_ufovc,'file')
            delete(fpath_ufovc);
        end
        fprintf('Saving %s ...',fpath_ufovc);
        save(fpath_ufovc,'msldemc_imUFOVmask','msldemc_imUFOVhdr',...
            'msldemc_imUFOVxynn','mastcam_msldemc_nn_UFOVmask',...
            'basename_msldem','dirpath_msldem');
        fprintf('\nDone.\n');
    end
    
else
    load(fpath_ufovc,'msldemc_imUFOVmask','msldemc_imUFOVhdr',...
        'msldemc_imUFOVxynn','mastcam_msldemc_nn_UFOVmask',...
        'basename_msldem','dirpath_msldem');
    if ~strcmpi(basename_msldem,MSLDEMdata.basename)
        error('MSLDEM used for the cache is different from that of input.');
    end

end

