function [msldemc_imUFOVmask,msldemc_imUFOVhdr,msldemc_imUFOVxynn,mastcam_msldemc_nn_UFOVmask,...
    mapper_mastcam2msldemc,mapper_msldemc2mastcam_mat,mapcell_msldemc2mastcam]...
    = mastcam_crop_msldemc_imUFOVmask(MSTproj,varargin)

global msl_env_vars
dirpath_cache = msl_env_vars.dirpath_cache;
save_file = true;
force = false;
load_cache_ifexist = true;
cache_vr = ''; % {'v0','v1'}

coordsys = 'NorthEastNadir';

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
            case 'COORDINATE_SYSTEM'
                coordsys = varargin{i+1};
                
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
            || isempty(MSTproj.mastcam_msldemc_nn)
        error('Processing required preceding computation is loaded');
    end
    
    %% Detect the margin
    valid_lines   = find(any(MSTproj.msldemc_imUFOVmask',1));
    lrnge         = [valid_lines(1), valid_lines(end)];
    len_vl        = lrnge(2)-lrnge(1)+1;
    valid_samples = find(any(MSTproj.msldemc_imUFOVmask,1));
    srnge         = [valid_samples(1), valid_samples(end)];
    len_vs        = srnge(2)-srnge(1)+1;
    
    if all(valid_lines) && all(valid_samples) && ~isempty(MSTproj.msldemc_imUFOVhdr)
        msldemc_imUFOVhdr = MSTproj.msldemc_imUFOVhdr;
        % msldemc_northing = msldemc_imUFOVhdr.y;
        % msldemc_easting  = msldemc_imUFOVhdr.x;
    else   

        l1   = MSTproj.msldemc_imFOVhdr.line_offset+lrnge(1);
        % lend = MSTproj.msldemc_imFOVhdr.line_offset+lrnge(2);
        s1   = MSTproj.msldemc_imFOVhdr.sample_offset+srnge(1);
        % send = MSTproj.msldemc_imFOVhdr.sample_offset+srnge(2);
        % msldemc_northing = MSLDEMdata.hdr.y(l1:lend);
        % msldemc_easting  = MSLDEMdata.hdr.x(s1:send);

        % Crop the mask
        msldemc_imUFOVmask = MSTproj.msldemc_imUFOVmask(lrnge(1):lrnge(2),srnge(1):srnge(2));

        % Updated header information (including offset and the size of the mask
        % image
        msldemc_imUFOVhdr = [];
        msldemc_imUFOVhdr.lines   = len_vl;
        msldemc_imUFOVhdr.samples = len_vs;
        msldemc_imUFOVhdr.line_offset   = l1-1;
        msldemc_imUFOVhdr.sample_offset = s1-1;
        % msldemc_imUFOVhdr.y = msldemc_northing;
        % msldemc_imUFOVhdr.x = msldemc_easting;
    end

    %% Read imxy, and fill NaN for the pixels that are not observed.
    %---------------------------------------------------------------------%
    % Get cam_C_geo
    %---------------------------------------------------------------------%
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
                    ' a subclass of MSLGaleDEMMosaic_v3' ...
                    ]);
            end
            [cmmdl_geo] = transform_CAHVOR_MODEL_ROVERNAV2IAUMARS(cmmdl, ...
                 rover_nav_coord,'Mars_Shape','Sphere');
        case {'IAU_MARS_ELLIPSOID'}
            if ~isa(MSLDEMdata,'MSLGaleMosaicRadius_v3')
                error([ ...
                    'MSLDEMdata needs to be an object of MSLGaleMosaicRadius_v3,' ...
                    ' a subclass of MSLGaleDEMMosaic_v3' ...
                    ]);
            end
            [cmmdl_geo] = transform_CAHVOR_MODEL_ROVERNAV2IAUMARS(cmmdl, ...
                 rover_nav_coord,'Mars_Shape','Ellipsoid');
        otherwise
            error('Undefined COORDINATE_SYSTEM %s',coordsys);
    end
    cmmdl_geo.get_image_plane_unit_vectors();
    
    s1   = msldemc_imUFOVhdr.sample_offset + 1;
    send = msldemc_imUFOVhdr.sample_offset + msldemc_imUFOVhdr.samples;
    l1   = msldemc_imUFOVhdr.line_offset   + 1;
    lend = msldemc_imUFOVhdr.line_offset   + msldemc_imUFOVhdr.lines;
    
    switch upper(coordsys)
        case {'NEE','NORTHEASTELEVATION','NORTHEASTNADIR','NENADIR'}
            msldemc_northing = MSLDEMdata.northing(l1:lend);
            msldemc_easting  = MSLDEMdata.easting(s1:send);
            switch mastcamdata_obj.Linearization
                case 1
                    tic; [msldemc_imUFOVx,msldemc_imUFOVy] = cahv_get_imxy_MSLDEM_mex(...
                        MSLDEMdata.imgpath,MSLDEMdata.hdr,msldemc_imUFOVhdr,...
                        msldemc_northing,msldemc_easting,msldemc_imUFOVmask,cmmdl_geo); toc;
                    % msldemc_imUFOVxy = cat(3,dem_imx,dem_imy);
                case 0
                    tic; [msldemc_imUFOVx,msldemc_imUFOVy] = cahvor_get_imxy_msldem_mex(...
                        MSLDEMdata.imgpath,MSLDEMdata.hdr,msldemc_imUFOVhdr,...
                        msldemc_northing,msldemc_easting,msldemc_imUFOVmask,cmmdl_geo); toc;
                    % msldemc_imUFOVxy = cat(3,dem_imx,dem_imy);

                otherwise
                    error('Linearization %d is not supported',mastcamdata_obj.Linearization);
            end
        case {'IAU_MARS_SPHERE','IAU_MARS_ELLIPSOID'}
            msldemc_latitude  = deg2rad(MSLDEMdata.latitude(l1:lend));
            msldemc_longitude = deg2rad(MSLDEMdata.longitude(s1:send));
            switch mastcamdata_obj.Linearization
                case 1
                    tic; [msldemc_imUFOVx,msldemc_imUFOVy] = cahv_iaumars_get_imxy_MSLDEM_mex(...
                        MSLDEMdata.imgpath,MSLDEMdata.hdr, MSLDEMdata.OFFSET, ...
                        msldemc_imUFOVhdr,...
                        msldemc_latitude,msldemc_longitude,msldemc_imUFOVmask,cmmdl_geo); toc;
                case 0
                    tic; [msldemc_imUFOVx,msldemc_imUFOVy] = cahvor_iaumars_get_imxy_MSLDEM_mex(...
                        MSLDEMdata.imgpath,MSLDEMdata.hdr, MSLDEMdata.OFFSET, ...
                        msldemc_imUFOVhdr,...
                        msldemc_latitude,msldemc_longitude,msldemc_imUFOVmask,cmmdl_geo); toc;

                otherwise
                    error('Linearization %d is not supported',mastcamdata_obj.Linearization);
            end
            
        otherwise
            error('Undefined COORDINATE_SYSTEM %s',coordsys);
    end
    % clear dem_imx dem_imy;
    % msldemc_imUFOVxy = MSTproj.msldemc_imFOVxy(lrnge(1):lrnge(2),srnge(1):srnge(2),:);
    
    % mm = (msldemc_imUFOVmask==0);
    % for i=1:2
    %     msldemc_imUFOVtmp = msldemc_imUFOVxy(:,:,i);
    %     msldemc_imUFOVtmp(mm) = nan;
    %     msldemc_imUFOVxy(:,:,i) = msldemc_imUFOVtmp;
    % end
    
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

    % pixel correspondence of mastcam_msldemc_nn is updated for the cropped
    % image.
    mastcam_msldemc_nn_UFOVmask = MSTproj.mastcam_msldemc_nn;
    mastcam_msldemc_nn_UFOVmask(:,:,1) = mastcam_msldemc_nn_UFOVmask(:,:,1)-int32(srnge(1)-1);
    mastcam_msldemc_nn_UFOVmask(:,:,2) = mastcam_msldemc_nn_UFOVmask(:,:,2)-int32(lrnge(1)-1);
    
    %% Create mapper array
    mapper_mastcam2msldemc = msl_create_mapping_mastcam2msldemc_mex_v2(...
    msldemc_imUFOVxnn,msldemc_imUFOVynn,...
    mastcam_msldemc_nn_UFOVmask(:,:,1)-1,mastcam_msldemc_nn_UFOVmask(:,:,2)-1);

    [mapidx_msldemc2mastcam,mapcell_msldemc2mastcam]...
        = msl_create_mapping_msldemc2mastcam_mex_v2(...
    msldemc_imUFOVxnn,msldemc_imUFOVynn,...
    mastcam_msldemc_nn_UFOVmask(:,:,1)-1,mastcam_msldemc_nn_UFOVmask(:,:,2)-1);

    mapper_msldemc2mastcam_mat  = sparse(double(mapidx_msldemc2mastcam+1));

    msldemc_imUFOVxynn = cat(3,msldemc_imUFOVxnn,msldemc_imUFOVynn)+1;
    
    
    %% saving files
    if save_file
        basename_msldem = MSLDEMdata.basename;
        dirpath_msldem  = MSLDEMdata.dirpath;
        if exist(fpath_ufovc,'file')
            delete(fpath_ufovc);
        end
        fprintf('Saving %s ...',fpath_ufovc);
        save(fpath_ufovc,'msldemc_imUFOVmask','msldemc_imUFOVhdr',...
            'msldemc_imUFOVxynn','mastcam_msldemc_nn_UFOVmask',...
            'mapper_mastcam2msldemc','mapper_msldemc2mastcam_mat','mapcell_msldemc2mastcam',...
            'basename_msldem','dirpath_msldem');
        fprintf('\nDone.\n');
    end
    
else
    load(fpath_ufovc,'msldemc_imUFOVmask','msldemc_imUFOVhdr',...
        'msldemc_imUFOVxynn','mastcam_msldemc_nn_UFOVmask',...
        'mapper_mastcam2msldemc','mapper_msldemc2mastcam_mat','mapcell_msldemc2mastcam',...
        'basename_msldem','dirpath_msldem');
    if ~strcmpi(basename_msldem,MSLDEMdata.basename)
        error('MSLDEM used for the cache is different from that of input.');
    end

end

