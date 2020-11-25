function [hdr_info] = mastcam_extract_imghdr_from_lbl(lbl)
% [hdr_info] = mastcam_extract_imghdr_from_lbl(lbl_info)
%   extract header information (envi format) from CRISM LABEL file
%  Input Parameters
%   lbl: struct of LABEL file
%  Output Parameters
%   hdr_info: struct of header in envi format, if no image is found, [] is
%             returend.

[ obj_file_image ] = lbl.OBJECT_IMAGE;

if isempty(obj_file_image)
    hdr_info = [];
else
    hdr_info = [];
    hdr_info.samples = lbl.OBJECT_IMAGE.LINE_SAMPLES;
    hdr_info.lines = lbl.OBJECT_IMAGE.LINES;
    hdr_info.bands = lbl.OBJECT_IMAGE.BANDS;
    
    switch upper(lbl.OBJECT_IMAGE.SAMPLE_TYPE)
        case 'PC_REAL'
            hdr_info.data_type = 4;
            hdr_info.byte_order = 0;
        case 'MSB_UNSIGNED_INTEGER'
            hdr_info.byte_order = 1;
            switch lbl.OBJECT_IMAGE.SAMPLE_BITS
                case 16
                    hdr_info.data_type = 12;
                case 8
                    hdr_info.data_type = 1;
                otherwise
                    error('Undefined "lbl.OBJECT_IMAGE.SAMPLE_BITS"');
            end
        case 'UNSIGNED_INTEGER'
            hdr_info.byte_order = 1;
            hdr_info.data_type = 1;
        case 'MSB_INTEGER'
            hdr_info.byte_order = 1;
            switch lbl.OBJECT_IMAGE.SAMPLE_BITS
                case 16
                    hdr_info.data_type = 2;
                otherwise
                    error('Undefined "lbl.OBJECT_IMAGE.SAMPLE_BITS"');
            end
        otherwise
            error('The data type: %s is not supported.',lbl.OBJECT_IMAGE.SAMPLE_TYPE);
    end

    switch upper(lbl.OBJECT_IMAGE.BAND_STORAGE_TYPE)
        case 'LINE_INTERLEAVED'
            hdr_info.interleave = 'bil';
        case 'BAND_SEQUENTIAL'
            hdr_info.interleave = 'bsq';
        otherwise
            error('The interleave: %s is not supported.',lbl.OBJECT_IMAGE.BAND_STORAGE_TYPE);
    end
    
    hdr_info.header_offset = 0;
    % hdr_info.header_offset = img_obj.RECORD_BYTES;
    
    if isfield(lbl.OBJECT_IMAGE,'BAND_NAME')
        hdr_info.band_names = obj_file_image.OBJECT_IMAGE.BAND_NAME;
    end
    
    hdr_info.data_ignore_value = lbl.OBJECT_IMAGE.INVALID_CONSTANT;

end




end