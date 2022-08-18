function [dem_img] = msldem_lazyenvireadRect(MSLDEMdata,...
   sample_offset,line_offset,samplesc,linesc,varargin)
% [dem_img] = msldem_lazyenvireadRect(MSLDEMdata,...
%    sample_offset,line_offset,samplesc,linesc,varargin)
% INPUTS
%   MSLDEMdata: HSI class object
%   sample_offset, line_offset: samples and lines to be offset
%   samplesc, linesc: size of the rectangle
% OUTPUTS
%   dem_img: array [linesc x samplesc]

%    

precision = 'double';
mex_mode = 1;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'PRECISION'
                precision = varargin{i+1};
            case 'MEX_MODE'
                mex_mode = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end


switch mex_mode
    case 0
        MSLDEMdata.fopen_img();
        typeName = 'single'; machine = 'ieee-le'; size_type = 4;
        fseek(MSLDEMdata.fid_img,size_type*(MSLDEMdata.hdr.samples*line_offset),-1);
        skipbytes_l = size_type*sample_offset;
        skipbytes_r = size_type*(MSLDEMdata.hdr.samples-sample_offset-samplesc);
        dem_img = nan(linesc,samplesc,typeName);
        for li = 1:linesc 
            fseek(MSLDEMdata.fid_img,skipbytes_l,0);
            dem_img(li,:) = fread(MSLDEMdata.fid_img,samplesc,typeName,0,machine);
            fseek(MSLDEMdata.fid_img,skipbytes_r,0);
        end
        MSLDEMdata.fclose_img();
    case 1
        [dem_img] = msldem_lazyenvireadRect_mex(MSLDEMdata.imgpath,MSLDEMdata.hdr,...
            sample_offset,line_offset,samplesc,linesc);
        dem_img = dem_img';
end

switch precision
    case 'double'
        dem_img = double(dem_img);
    case 'single'
end

dem_img(dem_img==MSLDEMdata.hdr.data_ignore_value) = nan;

end
