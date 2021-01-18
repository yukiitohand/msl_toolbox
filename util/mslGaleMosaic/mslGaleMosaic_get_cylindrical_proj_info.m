function [proj_info] = mslGaleMosaic_get_cylindrical_proj_info(lbl)

map_scale =  lbl.OBJECT_IMAGE_MAP_PROJECTION.MAP_SCALE.value;

map_resol = lbl.OBJECT_IMAGE_MAP_PROJECTION.MAP_RESOLUTION.value;
S0 = lbl.OBJECT_IMAGE_MAP_PROJECTION.SAMPLE_PROJECTION_OFFSET.value;
first_pxl = lbl.OBJECT_IMAGE_MAP_PROJECTION.SAMPLE_FIRST_PIXEL;
% last_pxl  = lbl.OBJECT_IMAGE_MAP_PROJECTION.SAMPLE_LAST_PIXEL;
ctr_lon = lbl.OBJECT_IMAGE_MAP_PROJECTION.CENTER_LONGITUDE.value;
% longitudes = ((first_pxl:last_pxl) - S0-1) ./ map_resol + ctr_lon;
lon1 = (first_pxl - S0-1) ./ map_resol + ctr_lon;

L0 = lbl.OBJECT_IMAGE_MAP_PROJECTION.LINE_PROJECTION_OFFSET.value;
first_pxl = lbl.OBJECT_IMAGE_MAP_PROJECTION.LINE_FIRST_PIXEL;
% last_pxl  = lbl.OBJECT_IMAGE_MAP_PROJECTION.LINE_LAST_PIXEL;
ctr_lat = lbl.OBJECT_IMAGE_MAP_PROJECTION.CENTER_LATITUDE.value;
% latitudes = (L0-(first_pxl:last_pxl)+1) ./ map_resol + ctr_lat;
lat1 = (L0-first_pxl+1) ./ map_resol + ctr_lat;


proj_info = MarsSphereEquiRectangularProj('Radius',3396190,...
    'STANDARD_PARALLEL',0,'CenterLongitude',[],...
    'Latitude_of_origin',0,'Longitude_of_origin',0);
proj_info.rdlat = map_resol;
proj_info.rdlon = map_resol;
proj_info.map_scale_x = map_scale;
proj_info.map_scale_y = map_scale;
proj_info.set_lat1(lat1);
proj_info.set_lon1(lon1);
proj_info.longitude_range = [lbl.OBJECT_IMAGE_MAP_PROJECTION.WESTERNMOST_LONGITUDE.value lbl.OBJECT_IMAGE_MAP_PROJECTION.EASTERNMOST_LONGITUDE.value];
proj_info.latitude_range = [lbl.OBJECT_IMAGE_MAP_PROJECTION.MAXIMUM_LATITUDE.value lbl.OBJECT_IMAGE_MAP_PROJECTION.MINIMUM_LATITUDE.value];

end