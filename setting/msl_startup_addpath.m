function [] = msl_startup_addpath()
%-------------------------------------------------------------------------%
% % Automatically find the path to toolboxes
fpath_self = mfilename('fullpath');
[dirpath_self,filename] = fileparts(fpath_self);

idx_sep = strfind(dirpath_self,'msl_toolbox/setting');
toolbox_root_dir = dirpath_self(1:idx_sep-1);


%-------------------------------------------------------------------------%
% name of the directory of each toolbox
base_toolbox_dirname          = 'base';
envi_toolbox_dirname          = 'envi';
pds3_toolbox_dirname          = 'pds3_toolbox';
crism_toolbox_dirname         = 'crism_toolbox';
msl_toolbox_dirname           = 'msl_toolbox';

%-------------------------------------------------------------------------%
% base toolbox
base_toolbox_dir = [toolbox_root_dir base_toolbox_dirname '/'];
addpath(base_toolbox_dir);
% joinPath in base toolbox will be used in the following. "base" toolbox
% need to be loaded first. base/joinPath.m automatically determine the
% presence of trailing slash, so you do not need to worry it.


% envi toolbox
envi_toolbox_dir = joinPath(toolbox_root_dir, envi_toolbox_dirname);
addpath( ...
    envi_toolbox_dir                                        , ...
    joinPath(envi_toolbox_dir,'v2/')                        , ...
    joinPath(envi_toolbox_dir,'v3/')                        , ...
    joinPath(envi_toolbox_dir,'v3/lazy_mex/')               , ...
    joinPath(envi_toolbox_dir,'v3/lazy_mex/wrapper/')         ...
);

cmp_arch = computer('arch');
switch cmp_arch
    case 'maci64'
        % For Mac computers
        addpath(joinPath(envi_toolbox_dir,'v3/lazy_mex/build/maci64/'));
    case 'glnxa64'
        % For Linux/Unix computers with x86-64 architechture, Sorry, no binary for Windows
        addpath(joinPath(envi_toolbox_dir,'v3/lazy_mex/build/glnxa64/'));
    case 'win64'
        fprintf('manually compile mex files\n');
    otherwise
        error('%s is not supported',cmp_arch);
end


% pds3_toolbox
pds3_toolbox_dir = joinPath(toolbox_root_dir, pds3_toolbox_dirname);
addpath( ...
    pds3_toolbox_dir                                , ...
    joinPath(pds3_toolbox_dir,'generic/')           , ...
    joinPath(pds3_toolbox_dir,'generic/connect/')   , ...
    joinPath(pds3_toolbox_dir,'generic/readwrite/') , ...
    joinPath(pds3_toolbox_dir,'generic/setting/')   , ...
    joinPath(pds3_toolbox_dir,'generic/util/')        ...
);

% crism_toolbox
crism_toolbox_dir = joinPath(toolbox_root_dir, crism_toolbox_dirname);
if ~exist(crism_toolbox_dir,'dir')
    fprintf(['crism_toolbox not installed. You can optionally install crism_toolbox at\n' ...
        '    github.com/yukiitohand/crism_toolbox/\n']);
else
    addpath( ...
        crism_toolbox_dir                                      , ...
        joinPath(crism_toolbox_dir,'base/')                    , ...
        joinPath(crism_toolbox_dir,'base/atf_util/')           , ...
        joinPath(crism_toolbox_dir,'base/basename_util/')      , ...
        joinPath(crism_toolbox_dir,'base/connect/')            , ...
        joinPath(crism_toolbox_dir,'base/folder_resolver/')    , ...
        joinPath(crism_toolbox_dir,'base/lbl_util/')           , ...
        joinPath(crism_toolbox_dir,'base/readwrite/')          , ...
        joinPath(crism_toolbox_dir,'core/')                    , ...
        joinPath(crism_toolbox_dir,'library/')                 , ...
        joinPath(crism_toolbox_dir,'library/base/')            , ...
        joinPath(crism_toolbox_dir,'library/conv/')            , ...
        joinPath(crism_toolbox_dir,'library/folder_resolver/') , ...
        joinPath(crism_toolbox_dir,'library/util/') , ...
        joinPath(crism_toolbox_dir,'map/')                     , ...
        joinPath(crism_toolbox_dir,'setting/')                 , ...
        joinPath(crism_toolbox_dir,'spice/')                   , ...
        joinPath(crism_toolbox_dir,'spice/cahv/')              , ...
        joinPath(crism_toolbox_dir,'spice/kernel/')            , ...
        joinPath(crism_toolbox_dir,'util/')                    , ...
        joinPath(crism_toolbox_dir,'util/ADRVS_util/')         , ...
        joinPath(crism_toolbox_dir,'util/BP_util/')            , ...
        joinPath(crism_toolbox_dir,'util/photocor/')             ...
    );
end

% msl_toolbox
msl_toolbox_dir = joinPath(toolbox_root_dir, msl_toolbox_dirname);
addpath( ...
    msl_toolbox_dir                                  , ...
    joinPath(msl_toolbox_dir, 'mastcam/')            , ...
    joinPath(msl_toolbox_dir, 'mastcam/base/')       , ...
    joinPath(msl_toolbox_dir, 'mastcam/core/')       , ...
    joinPath(msl_toolbox_dir, 'mastcam/INDEX/')      , ...
    joinPath(msl_toolbox_dir, 'mastcam/BANDPASS/')   , ...
    joinPath(msl_toolbox_dir, 'mastcam/projection/') , ...
    joinPath(msl_toolbox_dir, 'mastcam/projection/comp/')                , ...
    joinPath(msl_toolbox_dir, 'mastcam/projection/iau_mars/')            , ...
    joinPath(msl_toolbox_dir, 'mastcam/projection/iau_mars/comparison/') , ...
    joinPath(msl_toolbox_dir, 'mastcam/projection/mex/')        , ...
    joinPath(msl_toolbox_dir, 'mastcam/projection/mex/cahvor/') , ...
    joinPath(msl_toolbox_dir, 'mastcam/projection/mex/ufov/')   , ...
    joinPath(msl_toolbox_dir, 'mastcam/projection/polygons')    , ...
    joinPath(msl_toolbox_dir, 'mastcam/spectral/') , ...
    joinPath(msl_toolbox_dir, 'mastcam/util/')     , ...
    joinPath(msl_toolbox_dir, 'places/')      , ...
    joinPath(msl_toolbox_dir, 'places/maps/') , ...
    joinPath(msl_toolbox_dir, 'setting/')     , ...
    joinPath(msl_toolbox_dir, 'util/')        , ...
    joinPath(msl_toolbox_dir, 'util/mslGaleMosaic/')       , ...
    joinPath(msl_toolbox_dir, 'util/mslGaleMosaic/dem/')   , ...
    joinPath(msl_toolbox_dir, 'util/mslGaleMosaic/ortho/') , ...
    joinPath(msl_toolbox_dir, 'view/') ...
);

cmp_arch = computer('arch');
switch cmp_arch
    case 'maci64'
        % For Mac computers
        msl_mex_build_path = joinPath(msl_toolbox_dir,'mex_build/maci64/');
    case 'glnxa64'
        % For Linux/Unix computers with x86-64 architechture, Sorry, no support for Windows
        msl_mex_build_path = joinPath(msl_toolbox_dir,'mex_build/glnxa64/');
    case 'win64'
        msl_mex_build_path = joinPath(msl_toolbox_dir,'mex_build/win64/');
    otherwise
        error('%s is not supported',cmp_arch);
end

if exist(msl_mex_build_path,'dir')
    addpath(msl_mex_build_path);
else
    addpath(msl_mex_build_path);
    fprintf('Run msl_script_compile_all.m to compile C/MEX sources');
end
    
    

end