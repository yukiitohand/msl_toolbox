function [dirpath,subdir] = get_dirpath_mastcam(basename_mastcam)

global msl_env_vars
localrootDir = msl_env_vars.local_pds_msl_imaging_rootDir;
pds_msl_imaging_URL = msl_env_vars.pds_msl_imaging_URL;

mastcam_rootpath = joinPath(localrootDir,pds_msl_imaging_URL);
propMASTCAM = getProp_basenameMASTCAM(basename_mastcam);
[match_list] = mastcam_product_search(propMASTCAM);

subdir = joinPath(match_list(1).VOLUME_ID,match_list(1).PATH_NAME);
dirpath = joinPath(mastcam_rootpath,match_list(1).VOLUME_ID,match_list(1).PATH_NAME);

end
    