function [] = msl_startup_addpath()
%-------------------------------------------------------------------------%
% % Automatically find the path to toolboxes
fpath_self = mfilename('fullpath');
[dirpath_self,filename] = fileparts(fpath_self);

mtch = regexpi(dirpath_self,'(?<parent_dirpath>.*)/msl_toolbox[/]{0,1}','names');
toolbox_root_dir = mtch.parent_dirpath;


%-------------------------------------------------------------------------%
% name of the directory of each toolbox
base_toolbox_dirname          = 'base';
envi_toolbox_dirname          = 'envi';
pds3_toolbox_dirname          = 'pds3_toolbox';
crism_toolbox_dirname         = 'crism_toolbox';
msl_toolbox_dirname           = 'msl_toolbox';

%-------------------------------------------------------------------------%
pathCell = strsplit(path, pathsep);
%%
% base toolbox
base_toolbox_dir = [toolbox_root_dir '/' base_toolbox_dirname '/'];
addpath(base_toolbox_dir);
% joinPath in base toolbox will be used in the following. "base" toolbox
% need to be loaded first. base/joinPath.m automatically determine the
% presence of trailing slash, so you do not need to worry it.
if exist(base_toolbox_dir,'dir')
    if ~check_path_exist(base_toolbox_dir, pathCell)
        addpath(base_toolbox_dir);
    end
else
    warning([ ...
        'base toolbox is not detected. Download from' '\n' ...
        '   https://github.com/yukiitohand/base/'
        ]);
end

%%
% envi toolbox
envi_toolbox_dir = joinPath(toolbox_root_dir, envi_toolbox_dirname);

if exist(envi_toolbox_dir,'dir')
    if ~check_path_exist(envi_toolbox_dir, pathCell)
        run(joinPath(envi_toolbox_dir,'envi_startup_addpath'));
    end
else
    warning([ ...
        'envi toolbox is not detected. Download from' '\n' ...
        '   https://github.com/yukiitohand/envi/'
        ]);
end

%%
pds3_toolbox_dir = joinPath(toolbox_root_dir, pds3_toolbox_dirname);

if exist(pds3_toolbox_dir,'dir')
    if ~check_path_exist(pds3_toolbox_dir, pathCell)
        run(joinPath(pds3_toolbox_dir,'pds3_startup_addpath'));
    end
else
    warning([ ...
        'pds3_toolbox is not detected. Download from' '\n' ...
        '   https://github.com/yukiitohand/pds3_toolbox/'
        ]);
end

%%
% crism_toolbox
crism_toolbox_dir = joinPath(toolbox_root_dir, crism_toolbox_dirname);
if exist(crism_toolbox_dir,'dir')
    if ~check_path_exist(crism_toolbox_dir, pathCell)
        run(joinPath(crism_toolbox_dir,'crism_startup_addpath'));
    end
else
    warning([ ...
        'crism_toolbox is not detected. Download from' '\n' ...
        '   https://github.com/yukiitohand/crism_toolbox/'
        ]);
end

%% msl_toolbox
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
        % For Linux/Unix computers with x86-64 architechture
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
    fprintf('Run setting/msl_script_compile_all.m to compile C/MEX sources.\n');
end
    
    

end

%%
function [onPath] = check_path_exist(dirpath, pathCell)
    % pathCell = strsplit(path, pathsep, 'split');
    if ispc || ismac 
      onPath = any(strcmpi(dirpath, pathCell));
    else
      onPath = any(strcmp(dirpath, pathCell));
    end
end