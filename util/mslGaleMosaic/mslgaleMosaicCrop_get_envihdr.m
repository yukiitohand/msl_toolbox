function [hdr_c] = mslgaleMosaicCrop_get_envihdr(MSLGaledata,mslgalemosaiccrop_hdr,varargin)
% [hdr_c] = mslgaleMosaicCrop_get_envihdr(MSLGaledata,mslgalemosaiccrop_hdr,varargin)
%  INPUTS
%   MSLGaledata: MSLGaleMosaic_v3 obj, base MSLGaleMosaic data for which
%     cropping is performed.
%   mslgalemosaiccrop_hdr: struct having the four fields relating to
%   cropping information
%     samples, lines, sample_offset, line_offset
%  Required Parameters
%   'DATA_TYPE': data_type defined in ENVI HEADER
%   'BANDS'    : number of bands
%   'BAND_NAMES': band names field in ENVI HEADER
%  Optional Parameters
%   'DATA_IGNORE_VALUE': 
%  OUTPUTS
%   hdr_c: struct, envi header

%% Parsing INPUT
data_type         = [];
bands             = [];
band_names        = [];
data_ignore_value = [];
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'DATA_TYPE'
                data_type = varargin{i+1};
            case 'BANDS'
                bands = varargin{i+1};
            case 'BAND_NAMES'
                band_names = varargin{i+1};
            case 'DATA_IGNORE_VALUE'
                data_ignore_value = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

if isempty(bands)
    error('Provide Parameter "BANDS"');
end
if isempty(data_type)
    error('Provide Parameter "DATA_TYPE"');
end
if ~isempty(band_names)
    if ~iscell(band_names)
        error('Parameter "BAND_NAMES" needs to be a cell array');
    end
    if bands~=length(band_names)
        error('Parameter "BAND_NAMES" needs to be same length as bands');
    end
end

if ~isa(MSLGaledata,'MSLGaleMosaic_v3')
    error('MSLGaleDEMMosaic_v3 needs to be an object of MSLGaleDEMMosaic_v3');
end

if ~isfield(mslgalemosaiccrop_hdr,'samples')           ...
    || ~isfield(mslgalemosaiccrop_hdr,'lines')         ...
    || ~isfield(mslgalemosaiccrop_hdr,'sample_offset') ...
    || ~isfield(mslgalemosaiccrop_hdr,'line_offset')
        error([ ...
            'mslgalemosaiccrop_hdr is not right. Needs to have four fields:' ...
            '     samples'       ...
            '     lines'         ...
            '     sample_offset' ...
            '     line_offset'   ...
            ]);
end

%%
hdr_c = MSLGaledata.hdr;
hdr_c.samples = mslgalemosaiccrop_hdr.samples;
hdr_c.lines   = mslgalemosaiccrop_hdr.lines;
hdr_c.bands   = bands;
hdr_c.data_type = data_type;
if isempty(data_ignore_value)
    if isfield(hdr_c,'data_ignore_value')
        hdr_c = rmfield(hdr_c,'data_ignore_value');
    end
else
    hdr_c.data_ignore_value = data_ignore_value;
end
if isempty(band_names)
    if isfield(hdr_c,'band_names')
        hdr_c = rmfield(hdr_c,'band_names');
    end
else
    hdr_c.band_names = band_names;
end

hdr_c.map_info.mapx = MSLGaledata.easting(double(mslgalemosaiccrop_hdr.sample_offset+1));
hdr_c.map_info.mapy = MSLGaledata.northing(double(mslgalemosaiccrop_hdr.line_offset+1));

% offset parameters from the grid of the original MSLDEM Mosaic image
if isa(MSLGaledata,'MSLGaleDEMMosaic_v3')
    hdr_c.msldem_basename = MSLGaledata.basename;
    hdr_c.msldem_dirpath  = MSLGaledata.dirpath;
    hdr_c.msldemc_sample_offset = mslgalemosaiccrop_hdr.sample_offset;
    hdr_c.msldemc_line_offset   = mslgalemosaiccrop_hdr.line_offset;
elseif isa(MSLGaledata,'MSLGaleOrthoMosaic_v3')
    hdr_c.mslortho_basename = MSLGaledata.basename;
    hdr_c.mslortho_dirpath  = MSLGaledata.dirpath;
    hdr_c.mslorthoc_sample_offset = mslgalemosaiccrop_hdr.sample_offset;
    hdr_c.mslorthoc_line_offset   = mslgalemosaiccrop_hdr.line_offset;
else
    error('MSLGaledata is not appropriate');
end

end
