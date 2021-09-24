function [proj_info_from_base] = msldemc_get_proj_info_from_base(MSLDEMdata,raster_msldemc)
% [proj_info_from_base] = msldemc_get_proj_info_from_base(MSLDEMdata,raster_msldemc)
%   Get cylindical projection info class of ENVIRaster_MSLDEMC from its
%   base map
%  INPUTS
%   MSLDEMdata: object of MSLGaleDEMMosaic_v3 or MSLGaleMosaicRadius_v3
%   raster_msldemc : object of ENVIRasterMultBandMSLDEMCProj or 
%        ENVIRasterSingleLayerMSLDEMCProj
%  OUTPUTS
%   proj_info_from_base: SphereEquiRectangularProj object
% 

% Input check
validateattributes(MSLDEMdata,{'MSLGaleDEMMosaic_v3'},{},mfilename,'MSLDEMdata');
validateattributes(raster_msldemc, ...
    {'ENVIRasterMultBandMSLDEMCProj','ENVIRasterSingleLayerMSLDEMCProj'}, ...
    {},mfilename,'raster_msldemc');

proj_info_from_base = SphereEquiRectangularProj( ...
    'Radius'             , MSLDEMdata.proj_info.radius                , ...
    'STANDARD_PARALLEL'  , MSLDEMdata.proj_info.standard_parallel     , ...
    'CenterLongitude'    , MSLDEMdata.proj_info.center_longitude      , ...
    'Latitude_of_origin' , MSLDEMdata.proj_info.latitude_of_origin    , ...
    'Longitude_of_origin', MSLDEMdata.proj_info.longitude_of_origin     ...
    );
proj_info_from_base.rdlat = MSLDEMdata.proj_info.rdlat;
proj_info_from_base.rdlon = MSLDEMdata.proj_info.rdlon;
proj_info_from_base.map_scale_x = MSLDEMdata.proj_info.map_scale_x;
proj_info_from_base.map_scale_y = MSLDEMdata.proj_info.map_scale_y;
proj_info_from_base.set_lon1(MSLDEMdata.longitude(raster_msldemc.get_x_base(1)));
proj_info_from_base.set_lat1(MSLDEMdata.latitude(raster_msldemc.get_y_base(1)));
proj_info_from_base.longitude_range = MSLDEMdata.proj_info.get_lon_range(raster_msldemc.get_xrange_base());
proj_info_from_base.latitude_range  = MSLDEMdata.proj_info.get_lat_range(raster_msldemc.get_yrange_base());

end