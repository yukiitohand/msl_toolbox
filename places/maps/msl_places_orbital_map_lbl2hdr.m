function [hdr_info] = msl_places_orbital_map_lbl2hdr(lbl_info)
% [hdr_info] = msl_places_orbital_map_lbl2hdr(lbl_info)
%   extract header information (envi format) from MSL orbital_dem lbl
%  Input Parameters
%   lbl: struct of LABEL file
%  Output Parameters
%   hdr_info: struct of header in envi format, if no image is found, [] is
%             returend.



hdr_info = [];
hdr_info.samples = lbl_info.OBJECT_IMAGE.LINE_SAMPLES;
hdr_info.lines   = lbl_info.OBJECT_IMAGE.LINES;
hdr_info.bands   = lbl_info.OBJECT_IMAGE.BANDS;

[hdr_info.data_type,hdr_info.byte_order] = pds3_stsb2envihdr_dtbo(...
    lbl_info.OBJECT_IMAGE.SAMPLE_TYPE,lbl_info.OBJECT_IMAGE.SAMPLE_BITS);

[hdr_info.interleave] = pds3_bst2envihdr_interleave(...
    lbl_info.OBJECT_IMAGE.BAND_STORAGE_TYPE);

hdr_info.header_offset = lbl_info.RECORD_BYTES * (lbl_info.POINTER_IMAGE{2}-1);

if isfield(lbl_info.OBJECT_IMAGE,'BAND_NAME')
    hdr_info.band_names = obj_file_image.OBJECT_IMAGE.BAND_NAME;
end

% next, loading mapping information
map_info = [];
map_info.projection = lbl_info.GROUP_SURFACE_PROJECTION_PARMS.MAP_PROJECTION_TYPE;
map_info.image_coords = [0.5 0.5];
map_info.mapx = lbl_info.GROUP_SURFACE_PROJECTION_PARMS.Y_AXIS_MINIMUM;
map_info.mapy = lbl_info.GROUP_SURFACE_PROJECTION_PARMS.X_AXIS_MAXIMUM;
map_info.dx  = lbl_info.GROUP_SURFACE_PROJECTION_PARMS.MAP_SCALE{1}.value;
map_info.dy  = lbl_info.GROUP_SURFACE_PROJECTION_PARMS.MAP_SCALE{2}.value;
map_info.datum = 'D_Mars_2000_Sphere';
map_info.units = 'Meters';

xi = 1:hdr_info.samples;
yi = 1:hdr_info.lines;

x = map_info.dx*(xi-map_info.image_coords(1))+map_info.mapx;
y = -map_info.dy*(yi-map_info.image_coords(2))+map_info.mapy;

hdr_info.x = x;
hdr_info.y = y;

end