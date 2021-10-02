% msl_script_compile_all.m
% Script for compiling all the necessary C/MEX codes
%
% Automatically find the path to msl_toolbox and pds3_toolbox
% If you want to use the separate installation of pds3_toolbox, manually
% edit "pds3_toolbox_path"

fpath_self = mfilename('fullpath');
[dirpath_self,filename] = fileparts(fpath_self);

idx_sep = strfind(dirpath_self,'msl_toolbox/setting');
dirpath_toolbox = dirpath_self(1:idx_sep-1);

pds3_toolbox_path = joinPath(dirpath_toolbox, 'pds3_toolbox/');
msl_toolbox_path  = joinPath(dirpath_toolbox,  'msl_toolbox/');


%% Prior checking if necessary files are accessible.
if ~exist(pds3_toolbox_path,'dir')
    error('pds3_toolbox does not exist. Get at github.com/yukiitohand/pds3_toolbox/');
end

if ~exist(msl_toolbox_path,'dir')
    error('msl_toolbox does not exist.');
end

% Sorry, these functions are currently not supported for MATLAB 9.3 or
% earlier. These requires another 
if verLessThan('matlab','9.4')
    error('You need MATLAB version 9.4 or newer');
end

%% Set source code paths and the output directory path.
pds3_toolbox_mex_include_path = joinPath(pds3_toolbox_path, 'mex_include');
msl_toolbox_mex_include_path  = joinPath(msl_toolbox_path,  'mex_include');

% minimal set of mex files for the computation of mastcam projection using
% CAHVOR model and IAU_Mars (spherical/ellipsoidal) coordinate system
source_filepaths_cahvor_iau_mars = { ...
    joinPath(msl_toolbox_path,'mastcam/projection/iau_mars/','cahvor_iaumars_get_imFOVmask_MSLDEM_scf2_mex.c')     , ...
    joinPath(msl_toolbox_path,'mastcam/projection/iau_mars/','cahvor_iaumars_proj_mastcam2MSLDEM_v6_mex.c')        , ...
    joinPath(msl_toolbox_path,'mastcam/projection/iau_mars/','iaumars_get_msldemtUFOVmask_ctr_L2PBK_LL0_M3_mex.c') , ...
    joinPath(msl_toolbox_path,'mastcam/projection/iau_mars/','iaumars_get_msldemtUFOVmask_d_L2PBK_LL0_M3_mex.c')   , ...
    joinPath(msl_toolbox_path,'mastcam/projection/iau_mars/','cahvor_iaumars_get_imxy_MSLDEM_mex.c')               , ...
    joinPath(msl_toolbox_path,'mastcam/projection/mex/ufov/','get_msldemcUFOVmask_wVrtxPxlmsldemtUFOVmask.c')
};

% CAHV model and IAU_Mars (spherical/ellipsoidal) coordinate system
% future implementation
% source_filepaths_cahv_iau_mars = { };

% CAHVOR model and equirectangular projective coordinate system
% future implementation
% source_filepaths_cahvor_nenadir = { };

% CAHV model and equirectangular projective coordinate system
% future implementation
% source_filepaths_cahv_nenadir = { };

% MASTCAM <-> MSLDEM mapper utilities
source_filepaths_mapper = { ...
    joinPath(msl_toolbox_path,'util/','msl_create_mapping_mastcam2msldemc_mex_v2.c') , ...
    joinPath(msl_toolbox_path,'util/','msl_create_mapping_msldemc2mastcam_mex_v2.c')   ...
};

source_filepaths = [source_filepaths_cahvor_iau_mars source_filepaths_mapper];

switch computer
    case 'MACI64'
        out_dir = joinPath(msl_toolbox_path,'mex_build','./maci64/');
    case 'GLNXA64'
        out_dir = joinPath(msl_toolbox_path,'mex_build','./glnxa64/');
    case 'PCWIN64'
        out_dir = joinPath(msl_toolbox_path,'mex_build','./pcwin64/');
    otherwise
        error('Undefined computer type %s.\n',computer);
end

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

%% Compile files one by one
for i=1:length(source_filepaths)
    filepaths = source_filepaths{i};
    fprintf('Compiling %s ...\n',filepath);
    mex(filepath, '-R2018a', ['-I' pds3_toolbox_mex_include_path], ...
        ['-I' msl_toolbox_mex_include_path]);
end

% If you manually compile mex codes, include the two directories 
%     pds3_toolbox_mex_include_path
%     msl_toolbox_mex_include_path
% using -I option.
% Also do not forget to add '-R2018a'

