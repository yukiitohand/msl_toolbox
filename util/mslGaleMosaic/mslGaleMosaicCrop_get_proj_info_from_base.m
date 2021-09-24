function [proj_info_from_base] = mslGaleMosaicCrop_get_proj_info_from_base(MSLGaledata,raster_mslgalec)
% [proj_info_from_base] = mslGaleMosaicCrop_get_proj_info_from_base(MSLGaledata,raster_mslgalec)
%   Get cylindical projection info class of ENVIRaster_MSLDEMC from its
%   base map
%  INPUTS
%   MSLGaledata: object of MSLGaleMosaic_v3
%   raster_mslgalec : object of ENVIRasterMultBandMSLDEMCProj or 
%        ENVIRasterSingleLayerMSLDEMCProj or ENVIRasterSingleLayerMSLOrthoCProj
%  OUTPUTS
%   proj_info_from_base: SphereEquiRectangularProj object
% 
ENVIRasterSingleLayerMSLOrthoCProj
% Input check
validateattributes(MSLGaledata,{'MSLGaleMosaic_v3'},{},mfilename,'MSLGaledata');
validateattributes(raster_mslgalec, ...
    {'ENVIRasterMultBandMSLDEMCProj','ENVIRasterSingleLayerMSLDEMCProj', ...
     'ENVIRasterSingleLayerMSLOrthoCProj'}, ...
    {},mfilename,'raster_msldemc');

proj_info_from_base = SphereEquiRectangularProj( ...
    'Radius'             , MSLGaledata.proj_info.radius                , ...
    'STANDARD_PARALLEL'  , MSLGaledata.proj_info.standard_parallel     , ...
    'CenterLongitude'    , MSLGaledata.proj_info.center_longitude      , ...
    'Latitude_of_origin' , MSLGaledata.proj_info.latitude_of_origin    , ...
    'Longitude_of_origin', MSLGaledata.proj_info.longitude_of_origin     ...
    );
proj_info_from_base.rdlat = MSLGaledata.proj_info.rdlat;
proj_info_from_base.rdlon = MSLGaledata.proj_info.rdlon;
proj_info_from_base.map_scale_x = MSLGaledata.proj_info.map_scale_x;
proj_info_from_base.map_scale_y = MSLGaledata.proj_info.map_scale_y;
proj_info_from_base.set_lon1(MSLGaledata.longitude(raster_mslgalec.get_x_base(1)));
proj_info_from_base.set_lat1(MSLGaledata.latitude(raster_mslgalec.get_y_base(1)));
proj_info_from_base.longitude_range = MSLGaledata.proj_info.get_lon_range(raster_mslgalec.get_xrange_base());
proj_info_from_base.latitude_range  = MSLGaledata.proj_info.get_lat_range(raster_mslgalec.get_yrange_base());

end