function [msldemc_imUFOVmask] = mastcam_getUFOVwMSLDEM_mexw(MSTproj,varargin)
%[msldemc_imUFOVmask] = mastcam_getUFOVwMSLDEM_mexw(MSTproj,varargin)
%   Get image Unobstructed Field of View (UFOV) with the precomputed 
%   imageFOV mask. 
%  INPUTS
%   MSTproj: object of MASTCAMCameraProjectionMSLDEM
%  OUTPUTS
%   msldemc_imUFOVmask: int8 2 dimensonal array. The size is defined in 
%   MSTproj.msldemc_imFOVhdr.
%  OPTIONAL PARAMETERS
%   "BORDER_ASSESS_OPT: string, combination of 'dcl'
%     (default) 'd'
%   "DIRPATH_CACHE": string, 
%     (default) msl_env_vars.dirpath_cache
%   "UPDATE_CACHE": boolean, force updating cache files or not
%     (default) false

global msl_env_vars
dirpath_cache = msl_env_vars.dirpath_cache;

border_assess_opt = 'd';
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
            case 'BORDER_ASSESS_OPT'
                border_assess_opt = lower(varargin{i+1});
            
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
%-------------------------------------------------------------------------%
% Get the size of the mastcam image
%-------------------------------------------------------------------------%
L_im = mastcamdata_obj.L_im; S_im = mastcamdata_obj.S_im;
% switch class(mastcamdata_obj)
%     case 'MASTCAMdata'
%         L_im = mastcamdata_obj.hdr.lines; S_im = mastcamdata_obj.samples;
%     case 'MASTCAMgroup_eye'
%         L_im = mastcamdata_obj.L_im; S_im = mastcamdata_obj.S_im;
%     otherwise
%         error('second input needs to be either MASTCAMdata or MASTCAMgroup_eye class.');
% end

%-------------------------------------------------------------------------%
% Get cam_C_geo
%-------------------------------------------------------------------------%
cmmdl = mastcamdata_obj.CAM_MDL;
rover_nav_coord = mastcamdata_obj.ROVER_NAV;
cmmdl_geo = transform_CAHVOR_MODEL_wROVER_NAV(cmmdl,rover_nav_coord);
C_geo = cmmdl_geo.C;

%-------------------------------------------------------------------------%
% Create label name for cache files
%-------------------------------------------------------------------------%
[basename_cache] = mastcam_create_basename_cache(mastcamdata_obj);
msldemc_imFOVhdr = MSTproj.msldemc_imFOVhdr;


%%
%-------------------------------------------------------------------------%
% Evaluate the center of each pixel
%-------------------------------------------------------------------------%
basename_imUFOVmask_ctr = sprintf('%s_imUFOVmask_ctr_%s.mat',basename_cache,cache_vr);
fpath_imUFOVmask_ctr = joinPath(dirpath_cache,basename_imUFOVmask_ctr);
[flg_reproc] = doyouwanttoprocess(fpath_imUFOVmask_ctr,force,load_cache_ifexist);
if flg_reproc
    %---------------------------------------------------------------------%
    % Get DEM image
    %---------------------------------------------------------------------%

    tic; msldemc_img = msldem_lazyenvireadRect(MSLDEMdata,...
               MSTproj.msldemc_imFOVhdr.sample_offset,...
               MSTproj.msldemc_imFOVhdr.line_offset,...
               MSTproj.msldemc_imFOVhdr.samples,...
               MSTproj.msldemc_imFOVhdr.lines,...
               'precision','double'); toc;



    l1 = MSTproj.msldemc_imFOVhdr.line_offset+1;
    lend = MSTproj.msldemc_imFOVhdr.line_offset+MSTproj.msldemc_imFOVhdr.lines;
    s1 = MSTproj.msldemc_imFOVhdr.sample_offset+1;
    send = MSTproj.msldemc_imFOVhdr.sample_offset+MSTproj.msldemc_imFOVhdr.samples;
    msldemc_northing = MSLDEMdata.hdr.y(l1:lend);
    msldemc_easting  = MSLDEMdata.hdr.x(s1:send);

    msldemc_northing = msldemc_northing - C_geo(1);
    msldemc_easting = msldemc_easting - C_geo(2);
    msldemc_img = -msldemc_img-C_geo(3);
    
    % computation
    
    tic; [ msldemc_imUFOVmask_ctr ] = find_hidden_mastcamMSLDEM_v6_mex(...
        msldemc_img,...0
        msldemc_northing,...1
        msldemc_easting,...2
        MSTproj.msldemc_imFOVxy(:,:,1),...3
        MSTproj.msldemc_imFOVxy(:,:,2),...4
        MSTproj.msldemc_imFOVmask,...5
        S_im,L_im); toc;...6,7
    
    if save_file
        basename_msldem = MSLDEMdata.basename;
        dirpath_msldem  = MSLDEMdata.dirpath;
        if exist(fpath_imUFOVmask_ctr,'file')
            delete(fpath_imUFOVmask_ctr);
        end
        fprintf('Saving %s ...',fpath_imUFOVmask_ctr);
        save(fpath_imUFOVmask_ctr,'msldemc_imUFOVmask_ctr','msldemc_imFOVhdr',...
            'basename_msldem','dirpath_msldem');
        fprintf('\nDone.\n');
    end
else
    load(fpath_imUFOVmask_ctr,'msldemc_imUFOVmask_ctr','basename_msldem');
    if ~strcmpi(basename_msldem,MSLDEMdata.basename)
        error('MSLDEM used for the cache is different from that of input.');
    end
end    
msldemc_imUFOVmask = msldemc_imUFOVmask_ctr + int8(MSTproj.msldemc_imFOVmask_ctrnn);

%% Safeguarding
%-------------------------------------------------------------------------%
% Evaluate the vertexes of each pixel
%-------------------------------------------------------------------------%
if contains(border_assess_opt,'d')
    basename_imUFOVmask_wdedge = sprintf('%s_imUFOVmask_wdedge_%s.mat',basename_cache,cache_vr);
    fpath_imUFOVmask_wdedge = joinPath(dirpath_cache,basename_imUFOVmask_wdedge);
    [flg_reproc] = doyouwanttoprocess(fpath_imUFOVmask_wdedge,force,load_cache_ifexist);
    
    if flg_reproc
        tic; [dem_imxm,dem_imym] = get_imxyclm_MSLDEM_mex(...
            MSLDEMdata.imgpath,MSLDEMdata.hdr,MSTproj.msldemc_imFOVhdr,...
            msldemc_northing+C_geo(1),msldemc_easting+C_geo(2),...
            MSTproj.msldemc_imFOVmask,S_im,L_im,cmmdl_geo); toc;
        
        % Evaluate the corner of the pixels
        %  o---------o
        %  |  (c,l)  |
        %  |    X    |  
        %  |         |
        %  o---------o (c+1/2,l+1/2)

        tic; [ msldemc_imUFOVmask_d ] = find_hidden_edges_mastcamMSLDEM_v6_mex(...
            msldemc_img,...0
            msldemc_northing,...1
            msldemc_easting,...2
            MSTproj.msldemc_imFOVxy(:,:,1),...3
            MSTproj.msldemc_imFOVxy(:,:,2),...4
            MSTproj.msldemc_imFOVmask,...5
            S_im,L_im,...6.7
            dem_imxm,dem_imym); toc;

        clear dem_imxm dem_imym;

        % save msldemc_imUFOVmask_d.mat msldemc_imUFOVmask_d;

        tic; [ msldemc_imUFOVmask_wdedge ] = find_hidden_withEdges_mastcamMSLDEM(...
            msldemc_imUFOVmask_d,...0
            MSTproj.msldemc_imFOVmask); toc; ...1
            
        if save_file
            basename_msldem = MSLDEMdata.basename;
            dirpath_msldem  = MSLDEMdata.dirpath;
            if exist(fpath_imUFOVmask_wdedge,'file')
                delete(fpath_imUFOVmask_wdedge);
            end
            fprintf('Saving %s ...',fpath_imUFOVmask_wdedge);
            save(fpath_imUFOVmask_wdedge,'msldemc_imUFOVmask_wdedge','msldemc_imFOVhdr',...
                'basename_msldem','dirpath_msldem');
            fprintf('\nDone.\n');
        end
    else
        load(fpath_imUFOVmask_wdedge,'msldemc_imUFOVmask_wdedge','basename_msldem');
        
        if ~strcmpi(basename_msldem,MSLDEMdata.basename)
            error('MSLDEM used for the cache is different from that of input.');
        end
    end
    msldemc_imUFOVmask = msldemc_imUFOVmask_wdedge + msldemc_imUFOVmask;
end

%%
%-------------------------------------------------------------------------%
% Evaluate the center of column edges
%-------------------------------------------------------------------------%
if contains(border_assess_opt,'c')
    basename_imUFOVmask_wcedge = sprintf('%s_imUFOVmask_wcedge_%s.mat',basename_cache,cache_vr);
    fpath_imUFOVmask_wcedge = joinPath(dirpath_cache,basename_imUFOVmask_wcedge);
    [flg_reproc] = doyouwanttoprocess(fpath_imUFOVmask_wcedge,force,load_cache_ifexist);
    if flg_reproc
        
        % Evaluate the center of column edges
        %   ---------
        %  |  (c,l)  |
        %  o    X    o
        %  |         |
        %   ---------
        
        tic; [dem_imxm,dem_imym] = get_imxycm_MSLDEM_mex(...
            MSLDEMdata.imgpath,MSLDEMdata.hdr,MSTproj.msldemc_imFOVhdr,...
            msldemc_northing+C_geo(1),msldemc_easting+C_geo(2),...
            MSTproj.msldemc_imFOVmask,S_im,L_im,cmmdl_geo); toc;

        tic; [ msldemc_imUFOVmask_c ] = find_hidden_cedges_mastcamMSLDEM_v6_mex(...
            msldemc_img,...0
            msldemc_northing,...1
            msldemc_easting,...2
            MSTproj.msldemc_imFOVxy(:,:,1),...3
            MSTproj.msldemc_imFOVxy(:,:,2),...4
            MSTproj.msldemc_imFOVmask,...5
            S_im,L_im,...6.7
            dem_imxm,dem_imym); toc;

        clear dem_imxm dem_imym;

        % save msldemc_imUFOVmask_c.mat msldemc_imUFOVmask_c;

        tic; [ msldemc_imUFOVmask_wcedge ] = find_hidden_withcEdges_mastcamMSLDEM(...
            msldemc_imUFOVmask_c,...0
            MSTproj.msldemc_imFOVmask); toc; ...1
            
        if save_file
            basename_msldem = MSLDEMdata.basename;
            dirpath_msldem  = MSLDEMdata.dirpath;
            if exist(fpath_imUFOVmask_wcedge,'file')
                delete(fpath_imUFOVmask_wcedge);
            end
            fprintf('Saving %s ...',fpath_imUFOVmask_wcedge);
            save(fpath_imUFOVmask_wcedge,'msldemc_imUFOVmask_wcedge','msldemc_imFOVhdr',...
                'basename_msldem','dirpath_msldem');
            fprintf('\nDone.\n');
            
        end
    else
        load(fpath_imUFOVmask_wcedge,'msldemc_imUFOVmask_wcedge','basename_msldem');
        if ~strcmpi(basename_msldem,MSLDEMdata.basename)
            error('MSLDEM used for the cache is different from that of input.');
        end
    end
    msldemc_imUFOVmask = msldemc_imUFOVmask_wcedge + msldemc_imUFOVmask;
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
    basename_imUFOVmask_wledge = sprintf('%s_imUFOVmask_wledge_%s.mat',basename_cache,cache_vr);
    fpath_imUFOVmask_wledge = joinPath(dirpath_cache,basename_imUFOVmask_wledge);
    [flg_reproc] = doyouwanttoprocess(fpath_imUFOVmask_wledge,force,load_cache_ifexist);
    if flg_reproc
        
        tic; [dem_imxm,dem_imym] = get_imxylm_MSLDEM_mex(...
            MSLDEMdata.imgpath,MSLDEMdata.hdr,MSTproj.msldemc_imFOVhdr,...
            msldemc_northing+C_geo(1),msldemc_easting+C_geo(2),...
            MSTproj.msldemc_imFOVmask,S_im,L_im,cmmdl_geo); toc;


        tic; [ msldemc_imUFOVmask_l ] = find_hidden_ledges_mastcamMSLDEM_v6_mex(...
            msldemc_img,...0
            msldemc_northing,...1
            msldemc_easting,...2
            MSTproj.msldemc_imFOVxy(:,:,1),...3
            MSTproj.msldemc_imFOVxy(:,:,2),...4
            MSTproj.msldemc_imFOVmask,...5
            S_im,L_im,...6.7
            dem_imxm,dem_imym); toc;

        clear dem_imxm dem_imym;

        % save msldemc_imUFOVmask_l.mat msldemc_imUFOVmask_l;

        tic; [ msldemc_imUFOVmask_wledge ] = find_hidden_withlEdges_mastcamMSLDEM(...
            msldemc_imUFOVmask_l,...0
            MSTproj.msldemc_imFOVmask); toc; ...1
            
        if save_file
            basename_msldem = MSLDEMdata.basename;
            dirpath_msldem  = MSLDEMdata.dirpath;
            if exist(fpath_imUFOVmask_wledge,'file')
                delete(fpath_imUFOVmask_wledge);
            end
            fprintf('Saving %s ...',fpath_imUFOVmask_wledge);
            save(fpath_imUFOVmask_wledge,'msldemc_imUFOVmask_wledge','msldemc_imFOVhdr',...
                'basename_msldem','dirpath_msldem');
            fprintf('\nDone.\n');
        end
    else
        load(fpath_imUFOVmask_wledge,'msldemc_imUFOVmask_wledge','basename_msldem');
        if ~strcmpi(basename_msldem,MSLDEMdata.basename)
            error('MSLDEM used for the cache is different from that of input.');
        end
    end
    msldemc_imUFOVmask = msldemc_imUFOVmask_wledge + msldemc_imUFOVmask;
end

%% Summary
msldemc_imUFOVmask = int8(msldemc_imUFOVmask>0);

end

