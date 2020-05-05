function [hdr_info] = msl_orbital_dem_lbl2hdr(lbl_info)
% [hdr_info] = msl_orbital_dem_lbl2hdr(lbl_info)
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

switch upper(lbl_info.OBJECT_IMAGE.SAMPLE_TYPE)
    case 'PC_REAL'
        hdr_info.data_type = 4;
        hdr_info.byte_order = 0;
    case 'IEEE_REAL'
        hdr_info.data_type = 4;
        hdr_info.byte_order = 1;
    case 'MSB_UNSIGNED_INTEGER'
        hdr_info.byte_order = 1;
        switch obj_file_image.OBJECT_IMAGE.SAMPLE_BITS
            case 16
                hdr_info.data_type = 12;
            case 8
                hdr_info.data_type = 1;
            otherwise
                error('Undefined "OBJECT_IMAGE.SAMPLE_BITS"');
        end
    otherwise
        error('The data type: %s is not supported.',lbl_info.OBJECT_IMAGE.SAMPLE_TYPE);
end

hdr_info.header_offset = 0;
% hdr_info.header_offset = img_obj.RECORD_BYTES;

switch upper(lbl_info.OBJECT_IMAGE.BAND_STORAGE_TYPE)
    case 'LINE_INTERLEAVED'
        hdr_info.interleave = 'bil';
    case 'BAND_SEQUENTIAL'
        hdr_info.interleave = 'bsq';
    otherwise
        error('The interleave: %s is not supported.',lbl_info.OBJECT_IMAGE.BAND_STORAGE_TYPE);
end

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