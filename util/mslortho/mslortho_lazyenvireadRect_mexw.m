function [subimg] = mslortho_lazyenvireadRect_mexw(mslorthodata_obj,...
   sample_offset,line_offset,samplesc,linesc,varargin)
% [subimg] = mslortho_lazyenvireadRect_mexw(megdrdata_obj,...
%    sample_offset,line_offset,samplesc,linesc,varargin)
%   MEX wrapper.
% INPUTS
%   megdrdata_obj: MSLOrthoMosaicdata class object
%   sample_offset, line_offset: samples and lines to be offset
%   samplesc, linesc: size of the rectangle
% OUTPUTS
%   subimg: array [linesc x samplesc]
% s
% OPTIONAL PARAMETERS
%   PRECISION: char, string; data type of the output image.
%       'double', 'single', 'raw'
%      if 'raw', the data is returned with the original data type of the
%      image.
%      (default) 'double'
%     
%      
%    

precision = 'double';
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'PRECISION'
                precision = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end


[subimg] = mslortho_lazyenvireadRectUint8_mex(mslorthodata_obj.imgpath,mslorthodata_obj.hdr,...
            sample_offset,line_offset,samplesc,linesc);
subimg = subimg';

switch lower(precision)
    case 'double'
        subimg = double(subimg);
    case 'single'
        subimg = single(subimg);
    case 'uint8'
        subimg = uint8(subimg);
    case 'raw'
        
    otherwise
        error('Undefined precision %d',precision);
end

if isfield(mslorthodata_obj.hdr,'data_ignore_value')
    subimg(subimg==mslorthodata_obj.hdr.data_ignore_value) = nan;
end


end
