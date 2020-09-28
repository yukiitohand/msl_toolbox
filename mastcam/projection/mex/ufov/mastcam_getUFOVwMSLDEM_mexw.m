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
update_cache = false;
cache_vr = 'v1';
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'BORDER_ASSESS_OPT'
                border_assess_opt = lower(varargin{i+1});
            case 'DIRPATH_CACHE'
                dirpath_cache = varargin{i+1};
            case 'UPDATE_CACHE'
                update_cache = varargin{i+1};
            case 'CACHE_VERSION'
                cache_vr = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
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
if iscell(mastcamdata_obj.PRODUCT_ID)
    productID_repre = mastcamdata_obj.PRODUCT_ID{1};
else
    productID_repre = mastcamdata_obj.PRODUCT_ID;
end
propMASTCAMdata = getProp_basenameMASTCAM(productID_repre);
cam_code = propMASTCAMdata.cam_code;
if isnumeric(propMASTCAMdata.sol)
    sol = sprintf('%04d',propMASTCAMdata.sol);
end
if isnumeric(propMASTCAMdata.seq_id)
    seq_id = sprintf('%06d',propMASTCAMdata.seq_id);
end
site_id  = mastcamdata_obj.RMC.SITE;
drive_id = mastcamdata_obj.RMC.DRIVE;
pose_id  = mastcamdata_obj.RMC.POSE;
rsm_mc   = mastcamdata_obj.RMC.RSM;
basename_cache = sprintf('%s%s%s_SITE%03dDRIVE%04dPOSE%03dRSM%03d',sol,cam_code,seq_id,...
    site_id,drive_id,pose_id,rsm_mc);
msldemc_imFOVhdr = MSTproj.msldemc_imFOVhdr;

%-------------------------------------------------------------------------%
% Get DEM image
%-------------------------------------------------------------------------%

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
%%
%-------------------------------------------------------------------------%
% Evaluate the center of each pixel
%-------------------------------------------------------------------------%
basename_imUFOVmask_ctr = sprintf('%s_imUFOVmask_ctr_%s.mat',basename_cache,cache_vr);
fpath_imUFOVmask_ctr = joinPath(dirpath_cache,basename_imUFOVmask_ctr);
if ~update_cache && exist(fpath_imUFOVmask_ctr,'file')
    load(fpath_imUFOVmask_ctr,'msldemc_imUFOVmask_ctr');
else
    tic; [ msldemc_imUFOVmask_ctr ] = find_hidden_mastcamMSLDEM_v6_mex(...
        msldemc_img,...0
        msldemc_northing,...1
        msldemc_easting,...2
        MSTproj.msldemc_imFOVxy(:,:,1),...3
        MSTproj.msldemc_imFOVxy(:,:,2),...4
        MSTproj.msldemc_imFOVmask,...5
        S_im,L_im); toc;...6,7
    
    save(fpath_imUFOVmask_ctr,'msldemc_imUFOVmask_ctr','msldemc_imFOVhdr');
end    
msldemc_imUFOVmask = msldemc_imUFOVmask_ctr + int8(MSTproj.msldemc_imFOVmask_ctrnn);

%% Safeguarding
%-------------------------------------------------------------------------%
% Evaluate the vertexes of each pixel
%-------------------------------------------------------------------------%
if contains(border_assess_opt,'d')
    basename_imUFOVmask_wdedge = sprintf('%s_imUFOVmask_wdedge_%s.mat',basename_cache,cache_vr);
    % basename_imUFOVmask_wdedge = [basename_cache '_imUFOVmask_wdedge.mat'];
    fpath_imUFOVmask_wdedge = joinPath(dirpath_cache,basename_imUFOVmask_wdedge);
    if ~update_cache && exist(fpath_imUFOVmask_wdedge,'file')
        load(fpath_imUFOVmask_wdedge,'msldemc_imUFOVmask_wdedge');
    else
        tic; [dem_imxm,dem_imym] = get_imxyclm_MSLDEM_mex(...
            MSLDEMdata.imgpath,MSLDEMdata.hdr,MSTproj.msldemc_imFOVhdr,...
            msldemc_northing+C_geo(1),msldemc_easting+C_geo(2),...
            MSTproj.msldemc_imFOVmask,S_im,L_im,cmmdl_geo); toc;

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
        
        save(fpath_imUFOVmask_wdedge,'msldemc_imUFOVmask_wdedge','msldemc_imFOVhdr');
    
    end
    msldemc_imUFOVmask = msldemc_imUFOVmask_wdedge + msldemc_imUFOVmask;
end

%%
%-------------------------------------------------------------------------%
% Evaluate the center of column edges
%-------------------------------------------------------------------------%
if contains(border_assess_opt,'c')
    basename_imUFOVmask_wcedge = sprintf('%s_imUFOVmask_wcedge_%s.mat',basename_cache,cache_vr);
    % basename_imUFOVmask_wcedge = [basename_cache '_imUFOVmask_wcedge.mat'];
    fpath_imUFOVmask_wcedge = joinPath(dirpath_cache,basename_imUFOVmask_wcedge);
    if ~update_cache && exist(fpath_imUFOVmask_wcedge,'file')
        load(fpath_imUFOVmask_wcedge,'msldemc_imUFOVmask_wcedge');
    else
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
            
        save(fpath_imUFOVmask_wcedge,'msldemc_imUFOVmask_wcedge','msldemc_imFOVhdr');
    end
    msldemc_imUFOVmask = msldemc_imUFOVmask_wcedge + msldemc_imUFOVmask;
end
    
%%
%-------------------------------------------------------------------------%
% Evaluate the center of line edges
%-------------------------------------------------------------------------%
if contains(border_assess_opt,'l')
    basename_imUFOVmask_wledge = sprintf('%s_imUFOVmask_wledge_%s.mat',basename_cache,cache_vr);
    % basename_imUFOVmask_wledge = [basename_cache '_imUFOVmask_wledge.mat'];
    fpath_imUFOVmask_wledge = joinPath(dirpath_cache,basename_imUFOVmask_wledge);
    if ~update_cache && exist(fpath_imUFOVmask_wledge,'file')
        load(fpath_imUFOVmask_wledge,'msldemc_imUFOVmask_wledge');
    else
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
            
        save(fpath_imUFOVmask_wledge,'msldemc_imUFOVmask_wledge','msldemc_imFOVhdr');
    end
    msldemc_imUFOVmask = msldemc_imUFOVmask_wledge + msldemc_imUFOVmask;
end

msldemc_imUFOVmask = int8(msldemc_imUFOVmask>0);

end

