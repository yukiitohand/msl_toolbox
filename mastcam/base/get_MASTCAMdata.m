function [MASTCAMgroupList,RDRINDEX_matched2] = get_MASTCAMdata(varargin)


% USAGE:
% >> mascamdata_obj = get_MASTCAMdata('SOL',909,'CAM_CODE','M[R]{1}',...
%         'PRODUCT_TYPE','[DEC]{1}','DATA_PROC_CODE','DRLX','SEQ_ID',3977)

% optional variables: 
%  'BASENAME'
%  'SOL','CAM_CODE','PRODUCT_TYPE','DATA_PROC_CODE','SEQ_ID'
%  {'Download','Dwld'}, 'SITE_ID','DRIVE_ID','POSE_ID','RSM_ID'

global msl_env_vars
localrootDir = msl_env_vars.local_pds_msl_imaging_rootDir;
pds_msl_imaging_URL = msl_env_vars.pds_msl_imaging_URL;
mastcam_rootpath = joinPath(localrootDir,pds_msl_imaging_URL);

basename = {};

sol            = '';
cam_code       = '';
seq_id         = '';
command_num    = '';
cdpid_counter  = '';
unique_cdpid   = '';
product_type   = '';
gop            = '';
data_proc_code = '';
vr             = '';


site_id  = [];
drive_id = [];
pose_id  = [];
rsm_mc   = [];

dwld = 0;
vb   = 1;

rover_nav_ver = 'localized_interp';
rover_nav_mstcam_code = '';

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        if ~isempty(varargin{i+1})
            switch upper(varargin{i})
                case 'BASENAME'
                    basename = varargin{i+1};
                case 'SOL'
                    sol = varargin{i+1};
                case 'CAM_CODE'
                    cam_code = varargin{i+1};
                case 'SEQ_ID'
                    seq_id = varargin{i+1};
                case 'COMMAND_NUM'
                    command_num = varargin{i+1};
                case 'CDPID_COUNTER'
                    cdpid_counter = varargin{i+1};
                case 'UNIQUE_CDPID'
                    unique_cdpid = varargin{i+1};
                case 'PRODUCT_TYPE'
                    product_type = varargin{i+1};
                case 'GOP'
                    gop = varargin{i+1};
                case 'DATA_PROC_CODE'
                    data_proc_code = varargin{i+1};
                case 'VERSION'
                    vr = varargin{i+1};
                %-
                case 'SITE_ID'
                    site_id = varargin{i+1};
                case 'DRIVE_ID'
                    drive_id = varargin{i+1};
                case 'POSE_ID'
                    pose_id = varargin{i+1};
                case 'RSM_MC'
                    rsm_mc = varargin{i+1};
                case {'DOWNLOAD','DWLD'}
                    dwld = varargin{i+1};
                case 'VERBOSE'
                    vb = varargin{i+1};
                case {'ROVER_NAV_VERSION','ROVER_NAV_VER'}
                    rover_nav_ver = varargin{i+1};
                case 'ROVER_NAV_MSTCAM_CODE'
                    rover_nav_mstcam_code = varargin{i+1};
                otherwise
                    error('Unrecognized option: %s', varargin{i});   
            end
        end
    end
end

propMASTCAM_search = create_propMASTCAMbasename('SOL',sol,...
    'CAM_CODE',cam_code,'SEQ_ID',seq_id,'COMMAND_NUM',command_num,...
    'CDPID_COUNTER',cdpid_counter,'UNIQUE_CDPID',unique_cdpid,...
    'PRODUCT_TYPE',product_type,'GOP',gop,...
    'DATA_PROC_CODE',data_proc_code,'VERSION',vr);


if ~isempty(basename) && ischar(basename), basename = {basename}; end

if isempty(basename)
    [RDRINDEX_matched] = mastcam_product_search_v2(propMASTCAM_search);
    site_id_list  = cat(1,RDRINDEX_matched.ROVER_MOTION_COUNTER_SITE);
    drive_id_list = cat(1,RDRINDEX_matched.ROVER_MOTION_COUNTER_DRIVE);
    pose_id_list  = cat(1,RDRINDEX_matched.ROVER_MOTION_COUNTER_POSE);
    rsm_mc_list   = cat(1,RDRINDEX_matched.ROVER_MOTION_COUNTER_RSM);

    % match site_id
    if isempty(site_id)
        m_site_id = true(length(site_id_list),1);
    else
        m_site_id = (site_id_list==site_id);
    end
    
    % match drive_id
    if isempty(drive_id)
        m_drive_id = true(length(drive_id_list),1);
    else
        m_drive_id = (drive_id_list==drive_id);
    end
    
    % match pose_id
    if isempty(pose_id)
        m_pose_id = true(length(pose_id_list),1);
    else
        m_pose_id = (pose_id_list==pose_id);
    end
    
    % match remoste sensing mast motion counter
    if isempty(rsm_mc)
        m_rsm_mc = true(length(rsm_mc_list),1);
    else
        m_rsm_mc = (rsm_mc_list==rsm_mc);
    end
    
    idxBool_match = all([m_site_id,m_drive_id,m_pose_id,m_rsm_mc],2);
    
    RDRINDEX_matched2 = RDRINDEX_matched(idxBool_match);
    site_id_list_2  = site_id_list(idxBool_match);
    drive_id_list_2 = drive_id_list(idxBool_match);
    pose_id_list_2  = pose_id_list(idxBool_match);
    rsm_mc_list_2   = rsm_mc_list(idxBool_match);
    
    counter_list = [site_id_list_2 drive_id_list_2 pose_id_list_2 rsm_mc_list_2];
    [counter_unique,ia_cl,ic_cl] = unique(counter_list,'rows');
    if vb
        if isempty(counter_unique)
            fprintf('no image is selected.\n');
        elseif size(counter_unique,1)>1
            fprintf('Selected images\n');
            for i=1:size(counter_unique,1)
                match_idx = find(ic_cl==i);
                fprintf('(site, drive, pose, mc): PRODUCT ID\n');
                fprintf('(% 4d, % 5d, % 4d, %2d):\n', counter_unique(i,1),...
                    counter_unique(i,2),counter_unique(i,3),counter_unique(i,4));
                for im = match_idx
                    fprintf('%s\n',RDRINDEX_matched2(im).PRODUCT_ID);
                end
                fprintf('\n');
            end
        else % if only one combination is selected
            fprintf('(site, drive, pose, mc)\n');
            fprintf('(% 4d, % 5d, % 4d, %2d):\n', counter_unique(1),...
                    counter_unique(2),counter_unique(3),counter_unique(4));
            fprintf('Selected products:\n');
            for im = 1:length(RDRINDEX_matched2)
                fprintf('%s\n',RDRINDEX_matched2(im).PRODUCT_ID);
            end
        end
    end
    MASTCAMgroupList = [];
    for ic=1:size(counter_unique,1)
        match_idx = find(ic_cl==ic);
        MASTCAMgroup_u = MASTCAMgroup();
        for im=1:length(match_idx)
            i = match_idx(im);
            basename = RDRINDEX_matched2(i).PRODUCT_ID;
            subdir = joinPath(RDRINDEX_matched2(i).VOLUME_ID,RDRINDEX_matched2(i).PATH_NAME);
            dpath = joinPath(mastcam_rootpath,subdir);
            fnamelist = dir(dpath);
            if isempty(extractMatchedBasename_v2(basename,[{fnamelist.name}],'exact',0))
                pds_msl_imaging_downloader(subdir,'basenameptrn',basename,'dwld',dwld);
            end
            prop = getProp_basenameMASTCAM(basename);
            switch upper(prop.data_proc_code)
                case 'DRXX'
                    mastcamdata_i = MASTCAMdataDRXX(basename,dpath,...
                    'ROVER_NAV_VERSION',rover_nav_ver,'ROVER_NAV_MSTCAM_CODE',rover_nav_mstcam_code);
                case 'DRCX'
                    mastcamdata_i = MASTCAMdataDRCX(basename,dpath,...
                    'ROVER_NAV_VERSION',rover_nav_ver,'ROVER_NAV_MSTCAM_CODE',rover_nav_mstcam_code);
                case 'DRLX'
                    mastcamdata_i = MASTCAMdataDRLX(basename,dpath,...
                    'ROVER_NAV_VERSION',rover_nav_ver,'ROVER_NAV_MSTCAM_CODE',rover_nav_mstcam_code);
                case 'DRCL'
                    mastcamdata_i = MASTCAMdataDRCL(basename,dpath,...
                    'ROVER_NAV_VERSION',rover_nav_ver,'ROVER_NAV_MSTCAM_CODE',rover_nav_mstcam_code);
                otherwise
                    error('Undefined DATA_PROC_CODE %s',upper(prop.data_proc_code));
            end
                % mastcamdata_i = class_mstdata(basename,dpath,...
                %     'ROVER_NAV_VERSION',rover_nav_ver,'ROVER_NAV_MSTCAM_CODE',rover_nav_mstcam_code);
            MASTCAMgroup_u.append(mastcamdata_i);
        end
        MASTCAMgroupList = [MASTCAMgroupList,MASTCAMgroup_u];
    end
    
else
    error('This basename mode is not yet fully implemented');
    
end

end
