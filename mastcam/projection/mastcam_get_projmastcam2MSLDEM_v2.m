function [mstgrp_proj,msldemc_imFOVmask_pxlctrnn_obj,suppl_info] ...
    = mastcam_get_projmastcam2MSLDEM_v2(mstdata_obj,MSLDEMdata,msldemc_imFOVmask_obj,varargin)
% [mstgrp_proj,msldemc_imFOVmask_pxlctrnn_obj,suppl_info] ...
%     = mastcam_get_projmastcam2MSLDEM_v2(mstdata_obj,MSLDEMdata,msldemc_imFOVmask,varargin)
% Perform projection and get projection info of the center of the image
% pixels using pre-computed msldmc_imFOVmask
% INPUTS
%   mstdata_obj: either of 'MASTCAMdata' or 'MASTCAMgroup_wProcCode'
%   MSLDEMdata : MSLDEMGaleMosaic_v3
%   msldemc_imFOVmask_obj : ENVIRasterSingleLayerMSLDEMCProj
% OUTPUTS
%   mstgrp_obj : MASTCAMgroup_projection
%   msldemc_imFOVmask_pxlctrnn: ENVIRasterSingleLayerMSLDEMCProj
%   suppl_info : struct, supplemental information for the projection
% 
% OPTIONAL Parameters
%    'CACHE_DIRPATH': 
%       (default) msl_env_vars.dirpath_cache
%    'SAVE_FILE': 
%       (default) 1
%    'FORCE': boolean, whether to process forcefully
%       (default) 0
%    {'LOAD_CACHE_IFEXIST','LCIE'}
%       (default) 1
%    'VARARGIN_PROCESS'
%       cell array, Optional Parameters for 
%            proj_mastcam2MSLDEM_v5_mexw_v2
%       (default) {}
%
% 

global msl_env_vars
pdir_cache = msl_env_vars.dirpath_cache;
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
                pdir_cache = varargin{i+1};
                validateattributes(pdir_cache,{'char'},{},mfilename,'CACHE_DIRPATH');
            case 'SAVE_FILE'
                save_file = varargin{i+1};
                validateattributes(save_file,{'numeric','logical'},{'binary'},mfilename,'SAVE_FILE');
            case 'FORCE'
                force = varargin{i+1};
                validateattributes(force,{'numeric','logical'},{'binary'},mfilename,'FORCE');
            case {'LOAD_CACHE_IFEXIST','LCIE'}
                load_cache_ifexist = varargin{i+1};
                validateattributes(load_cache_ifexist,{'numeric','logical'}, ...
                    {'binary'},mfilename,'LOAD_CACHE_IFEXIST');
            case {'CACHE_VER','CACHE_VERSION'}
                cache_vr = varargin{i+1};
                validateattributes(cache_vr,{'char'},{},mfilename,'CACHE_VERSION');
                
            % ## PROCESSING OPTIONS #--------------------------------------
            case 'VARARGIN_PROCESS'
                varargin_proc = varargin{i+1};
                validateattributes(varargin_proc,{'cell'},{},mfilename,'VARARGIN_PROCESS');
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

if isempty(cache_vr)
   error('Please enter cache version with the option "CACHE_VER" or "CACHE_VERSION".'); 
end

% Input validation
validateattributes(mstdata_obj, ...
    {'MASTCAMdata','MASTCAMgroup_wProcCode'},{},mfilename,'mstdata_obj');
validateattributes(MSLDEMdata,{'MSLGaleDEMMosaic_v3'},{},mfilename,'MSLDEMdata');
validateattributes(msldemc_imFOVmask_obj, ...
    {'ENVIRasterSingleLayerMSLDEMCProj','ENVIRasterMultBandMSLDEMCProj'}, ...
    {},mfilename,'msldemc_imFOVmask_obj');

%%
[cache_dirname] = mastcam_get_cache_dirname(mstdata_obj);
dirpath_cache = joinPath(pdir_cache,cache_dirname);
if ~exist(dirpath_cache,'dir')
    mkdir(dirpath_cache);
end
[basename_cache_com] = mastcam_create_basename_cache(mstdata_obj);

% XYZ (Depending on the coordinate system)
% band 1: X, band 2: Y, band 3: Z
% double precision
basename_xyz = mastcam_get_basename_cache(basename_cache_com,'mastcam_XYZ',cache_vr);
[xyz_imgpath,xyz_hdrpath] = mastcam_get_cachefilepath(basename_xyz,dirpath_cache);

% Northing-Easting
% band 1: northing, band 2: easting
% double precision
basename_ne = mastcam_get_basename_cache(basename_cache_com,'mastcam_NE',cache_vr);
[ne_imgpath,ne_hdrpath] = mastcam_get_cachefilepath(basename_ne,dirpath_cache);

% Latitude-Longitude
% band 1: latitude [deg], band 2: longitude [deg]
% double precision
basename_latlon = mastcam_get_basename_cache(basename_cache_com,'mastcam_latlon',cache_vr);
[latlon_imgpath,latlon_hdrpath] = mastcam_get_cachefilepath(basename_latlon,dirpath_cache);

% Zenith (Depending on the coordinate system)
% band 1: zenith
% double precision
basename_zenith = mastcam_get_basename_cache(basename_cache_com,'mastcam_zenith',cache_vr);
[zenith_imgpath,zenith_hdrpath] = mastcam_get_cachefilepath(basename_zenith,dirpath_cache);

% Range of the center of pixels
% band 1: range
% double precision
basename_range = mastcam_get_basename_cache(basename_cache_com,'mastcam_range',cache_vr);
[range_imgpath,range_hdrpath] = mastcam_get_cachefilepath(basename_range,dirpath_cache);

% surface plane parameters ([n1 n2 n3]x = c)
% band 1: n1, band 2: n2, band 3: n3, band 4: c
% double precision
basename_surfplnc = mastcam_get_basename_cache(basename_cache_com,'mastcam_surfplnc',cache_vr);
[surfplnc_imgpath,surfplnc_hdrpath] = mastcam_get_cachefilepath(basename_surfplnc,dirpath_cache);

% emission angle
% band 1: emission angle
% double precision
basename_emiang = mastcam_get_basename_cache(basename_cache_com,'mastcam_emiang',cache_vr);
[emiang_imgpath,emiang_hdrpath] = mastcam_get_cachefilepath(basename_emiang,dirpath_cache);

% nearest neighbor pixel in msldem [x_msldem,y_msldem]
% band 1: x_msldem, band 2: y_msldem
% int32 precision
basename_nn_msldem = mastcam_get_basename_cache(basename_cache_com,'mastcam_nn_msldemc',cache_vr);
[nn_msldem_imgpath,nn_msldem_hdrpath] = mastcam_get_cachefilepath(basename_nn_msldem,dirpath_cache);

% reference triangle on which the intersection of the pixel center vector
% present.
% band 1: x_msldem, band 2: y_msldem, band 3: ti(triangle id {1,2})
% int32 precision
basename_ref_msldem = mastcam_get_basename_cache(basename_cache_com,'mastcam_ref_msldemc',cache_vr);
[ref_msldem_imgpath,ref_msldem_hdrpath] = mastcam_get_cachefilepath(basename_ref_msldem,dirpath_cache);

% msldemc_imFOVmask
% int8 type
basename_msldemc_imFOVmask_pxlctrnn = mastcam_get_basename_cache( ...
    basename_cache_com,'msldemc_imFOVmask_pxlctrnn',cache_vr);
[msldemc_imFOVmask_pxlctrnn_imgpath,msldemc_imFOVmask_pxlctrnn_hdrpath] ...
    = mastcam_get_cachefilepath(basename_msldemc_imFOVmask_pxlctrnn,dirpath_cache);

% supplement info
supplpath = joinPath(dirpath_cache,[basename_cache_com '_mastcam2MSLDEM_suppl.mat']);

% collect all 
cachepaths = {...
    xyz_imgpath,xyz_hdrpath,ne_imgpath,ne_hdrpath,                         ...
    latlon_imgpath,latlon_hdrpath,zenith_imgpath,zenith_hdrpath,           ...
    range_imgpath,range_hdrpath,surfplnc_imgpath,surfplnc_hdrpath,         ...
    emiang_imgpath,emiang_hdrpath,nn_msldem_imgpath,nn_msldem_hdrpath,   ...
    ref_msldem_imgpath,ref_msldem_hdrpath,                               ...
    msldemc_imFOVmask_pxlctrnn_imgpath,msldemc_imFOVmask_pxlctrnn_hdrpath, ...
    supplpath ...
    };

% Evaluate to perform the 
[flg_reproc] = doyouwanttoprocess(cachepaths,force,load_cache_ifexist);

%%
if flg_reproc
    [mastcam_XYZ,mastcam_latlon,mastcam_NE,mastcam_zenith, ...
        mastcam_ref_msldem,mastcam_range,mastcam_nn_msldem,...
        msldemc_imFOVmask_pxlctrnn,mastcam_emi,mastcam_surfplnc,option] ...
        = proj_mastcam2MSLDEM_v5_mexw_v2(mstdata_obj,MSLDEMdata,msldemc_imFOVmask_obj,varargin_proc{:});
    
    [hdr_xyz] = envi_create_envihdr_from_img(mastcam_XYZ, ...
        'description','{XYZ coordinate of the MASTCAM Image Pixel Centers. }',...
        'band_names',{'X','Y','Z'});
    hdr_xyz.coordsys = option.xyz_coordsys;
    hdr_xyz.coordsys_unit = 'meters';
    [hdr_ne]  = envi_create_envihdr_from_img(mastcam_NE, ...
        'description','{Northing and Easting of the MASTCAM Image Pixel Centers. }',...
        'band_names',{'northing','easting'});
    hdr_ne.coordsys_unit = 'meters';
    [hdr_latlon] = envi_create_envihdr_from_img(mastcam_latlon, ...
        'description','{Latitude and Longitude of the MASTCAM Image Pixel Centers. }',...
        'band_names',{'latitude [deg]','longitude [deg]'});
    [hdr_zenith] = envi_create_envihdr_from_img(mastcam_zenith);
    switch option.xyz_coordsys
        case 'IAU_MARS_SPHERE'
            if strcmpi(class(MSLDEMdata),'MSLGaleMosaicRadius_v3')
                hdr_zenith.zenith_type = 'Radius';
            else
                error('Unreasonable combination of the MSLDEM data class %s and coordinate system %s',...
                    class(MSLDEMdata),option.coordsys);
            end
        case 'NorthEastNadir'
            if strcmpi(class(MSLDEMdata),'MSLGaleMosaicRadius_v3')
                hdr_zenith.zenith_type = 'Radius';
                hdr_zenith.radius_offset = MSLDEMdata.OFFSET;
            elseif strcmpi(class(MSLDEMdata),'MSLGaleDEMMosaic_v3')
                hdr_zenith.zenith_type = 'Topography';
            else
                error('Unreasonable combination of the MSLDEM data class %s and coordinate system %s',...
                    class(MSLDEMdata),option.coordsys);
            end
        otherwise
            error('Undefined coordinate system %s',option.coordsys);
    end
    hdr_zenith.description = sprintf('{%s of the MASTCAM Image Pixel Centers. }',hdr_zenith.zenith_type);
    hdr_zenith.band_names  = {sprintf('%s',hdr_zenith.zenith_type)};
    hdr_zenith.zenith_unit = 'meters';
    
    [hdr_range]  = envi_create_envihdr_from_img(mastcam_range, ...
        'description','{Range of the MASTCAM Image Pixel Centers.[meters]}',...
        'band_names',{'Range [m]'});
    [hdr_emiang] = envi_create_envihdr_from_img(mastcam_emi, ...
        'description','{Emission angles [deg] of the MASTCAM Image Pixel Centers.[meters]}',...
        'band_names',{'Emission Angles [deg]'});
    [hdr_surfplnc]    = envi_create_envihdr_from_img(mastcam_surfplnc, ...
        'description',['{Surface plane parameters at the MASTCAM Image Pixel Centers. ', ...
        'First three bands are normal vectors, the last one is the plane constant.}'],...
        'band_names',{'n1','n2','n3','c'});
    hdr_surfplnc.coordsys = option.xyz_coordsys;
    
    [hdr_nn_msldem]  = envi_create_envihdr_from_img(mastcam_nn_msldem, ...
        'description', ...
        ['{Nearest Neighbor pixel index of the reference msldem cropped image', ...
         ' for mastcam image pixel centers. First index is 1 (MATLAB style).)}'], ...
        'band_names',{'sample','line'});
    hdr_nn_msldem.msldem_basename = msldemc_imFOVmask_obj.hdr.msldem_basename;
    hdr_nn_msldem.msldem_dirpath  = msldemc_imFOVmask_obj.hdr.msldem_dirpath;
    hdr_nn_msldem.index_style = 'MATLAB';
    
    [hdr_ref_msldem] = envi_create_envihdr_from_img(mastcam_ref_msldem, ...
        'description', ...
        ['{Index of the surface triangle for mastcam image pixel centers.' ...
         'First index is 1 (MATLAB-style).)}'], ...
         'band_names',{'sample index','line index','triangle index'});
    hdr_ref_msldem.msldem_basename = msldemc_imFOVmask_obj.hdr.msldem_basename;
    hdr_ref_msldem.msldem_dirpath  = msldemc_imFOVmask_obj.hdr.msldem_dirpath;
    hdr_ref_msldem.index_style = 'MATLAB';
    
    [data_type] = envihdr_get_data_type_from_precision('int8');
    hdr_msldemc_imFOVmask_pxlctrnn = mslgaleMosaicCrop_get_envihdr( ...
        MSLDEMdata,msldemc_imFOVmask_obj.chdr,'BANDS',1,'DATA_TYPE',data_type);
    
    if nargout>2
        suppl_info = struct('msldem_basename',MSLDEMdata.basename, ...
            'msldem_dirpath',MSLDEMdata.dirpath,'option',option);
    end
    

    if save_file
        envi_save_raster(mastcam_XYZ,hdr_xyz,xyz_imgpath,xyz_hdrpath);
        envi_save_raster(mastcam_NE,hdr_ne,ne_imgpath,ne_hdrpath);
        envi_save_raster(mastcam_latlon,hdr_latlon,latlon_imgpath,latlon_hdrpath);
        envi_save_raster(mastcam_zenith,hdr_zenith,zenith_imgpath,zenith_hdrpath);
        envi_save_raster(mastcam_range,hdr_range,range_imgpath,range_hdrpath);
        envi_save_raster(mastcam_surfplnc,hdr_surfplnc,surfplnc_imgpath,surfplnc_hdrpath);
        envi_save_raster(mastcam_emi,hdr_emiang,emiang_imgpath,emiang_hdrpath);
        envi_save_raster(mastcam_nn_msldem,hdr_nn_msldem,nn_msldem_imgpath,nn_msldem_hdrpath);
        envi_save_raster(mastcam_ref_msldem,hdr_ref_msldem,ref_msldem_imgpath,ref_msldem_hdrpath);
        envi_save_raster(msldemc_imFOVmask_pxlctrnn,hdr_msldemc_imFOVmask_pxlctrnn, ...
            msldemc_imFOVmask_pxlctrnn_imgpath,msldemc_imFOVmask_pxlctrnn_hdrpath);
        
        mstgrp_proj = MASTCAMgroup_projection();
        mstgrp_proj.XYZ = ENVIRasterMultBand(basename_xyz,dirpath_cache);
        mstgrp_proj.NE  = ENVIRasterMultBand(basename_ne,dirpath_cache);
        mstgrp_proj.latlon = ENVIRasterMultBand(basename_latlon,dirpath_cache);
        mstgrp_proj.zenith = ENVIRasterSingleLayer(basename_zenith,dirpath_cache);
        mstgrp_proj.range  = ENVIRasterSingleLayer(basename_range,dirpath_cache);
        mstgrp_proj.emiang = ENVIRasterSingleLayer(basename_emiang,dirpath_cache);
        mstgrp_proj.surfplnc   = ENVIRasterMultBand(basename_surfplnc,dirpath_cache);
        mstgrp_proj.nn_msldem  = ENVIRasterMultBand(basename_nn_msldem,dirpath_cache);
        mstgrp_proj.ref_msldem = ENVIRasterMultBand(basename_ref_msldem,dirpath_cache);
        msldemc_imFOVmask_pxlctrnn_obj = ENVIRasterSingleLayerMSLDEMCProj( ...
            basename_msldemc_imFOVmask_pxlctrnn,dirpath_cache);
        
        if exist(supplpath,'file'),delete(supplpath); end
        fprintf('Saving %s ...',supplpath);
        msldem_basename = MSLDEMdata.basename;
        msldem_dirpath  = MSLDEMdata.dirpath;
        save(supplpath,'msldem_basename','msldem_dirpath','option');
        fprintf('\nDone.\n');
        
    else
        mstgrp_proj = MASTCAMgroup_projection();
        mstgrp_proj.XYZ = ENVIRasterMultBand('','');
        mstgrp_proj.XYZ.hdr = hdr_xyz;
        mstgrp_proj.XYZ.img = mastcam_XYZ;
        mstgrp_proj.NE = ENVIRasterMultBand('','');
        mstgrp_proj.NE.hdr = hdr_ne;
        mstgrp_proj.NE.img = mastcam_NE;
        mstgrp_proj.latlon = ENVIRasterMultBand('','');
        mstgrp_proj.latlon.hdr = hdr_latlon;
        mstgrp_proj.latlon.img = mastcam_latlon;
        mstgrp_proj.zenith = ENVIRasterSingleLayer('','');
        mstgrp_proj.zenith.hdr = hdr_zenith;
        mstgrp_proj.zenith.img = mastcam_zenith;
        mstgrp_proj.range = ENVIRasterSingleLayer('','');
        mstgrp_proj.range.hdr = hdr_range;
        mstgrp_proj.range.img = mastcam_range;
        mstgrp_proj.emiang = ENVIRasterSingleLayer('','');
        mstgrp_proj.emiang.hdr = hdr_emiang;
        mstgrp_proj.emiang.img = mastcam_emi;
        mstgrp_proj.surfplnc = ENVIRasterMultBand('','');
        mstgrp_proj.surfplnc.hdr = hdr_surfplnc;
        mstgrp_proj.surfplnc.img = mastcam_surfplnc;
        mstgrp_proj.nn_msldem = ENVIRasterMultBand('','');
        mstgrp_proj.nn_msldem.hdr = hdr_nn_msldem;
        mstgrp_proj.nn_msldem.img = mastcam_nn_msldem;
        mstgrp_proj.ref_msldem = ENVIRasterMultBand('','');
        mstgrp_proj.ref_msldem.hdr = hdr_ref_msldem;
        mstgrp_proj.ref_msldem.img = mastcam_ref_msldem;
        
        msldemc_imFOVmask_pxlctrnn_obj = ENVIRasterSingleLayerMSLDEMCProj('','');
        msldemc_imFOVmask_pxlctrnn_obj.img  = msldemc_imFOVmask_pxlctrnn;
        msldemc_imFOVmask_pxlctrnn_obj.hdr  = hdr_msldemc_imFOVmask_pxlctrnn;
        msldemc_imFOVmask_pxlctrnn_obj.chdr = msldemc_imFOVmask_obj.chdr;
        msldemc_imFOVmask_pxlctrnn_obj.proj_info = envi_get_proj_info_SphereEquirect( ...
            msldemc_imFOVmask_pxlctrnn_obj.hdr);
        
    end
else
    
    mstgrp_proj = MASTCAMgroup_projection();
    mstgrp_proj.XYZ = ENVIRasterMultBand(basename_xyz,dirpath_cache);
    mstgrp_proj.NE  = ENVIRasterMultBand(basename_ne,dirpath_cache);
    mstgrp_proj.latlon = ENVIRasterMultBand(basename_latlon,dirpath_cache);
    mstgrp_proj.zenith = ENVIRasterSingleLayer(basename_zenith,dirpath_cache);
    mstgrp_proj.range  = ENVIRasterSingleLayer(basename_range,dirpath_cache);
    mstgrp_proj.emiang = ENVIRasterSingleLayer(basename_emiang,dirpath_cache);
    mstgrp_proj.surfplnc   = ENVIRasterMultBand(basename_surfplnc,dirpath_cache);
    mstgrp_proj.nn_msldem  = ENVIRasterMultBand(basename_nn_msldem,dirpath_cache);
    mstgrp_proj.ref_msldem = ENVIRasterMultBand(basename_ref_msldem,dirpath_cache);

    msldemc_imFOVmask_pxlctrnn_obj = ENVIRasterSingleLayerMSLDEMCProj( ...
            basename_msldemc_imFOVmask_pxlctrnn,dirpath_cache);
    
    
    suppl_info = load(supplpath);
    if ~strcmpi(suppl_info.msldem_basename,MSLDEMdata.basename)
        error('MSLDEM used for the cache is different from that of input.');
    end
            
end

end
