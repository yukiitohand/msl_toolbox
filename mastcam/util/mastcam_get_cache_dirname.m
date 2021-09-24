function [cache_dirname] = mastcam_get_cache_dirname(mastcamdata_obj,varargin)
% [cache_dirname] = mastcam_get_cache_dirname(mastcamdata_obj,varargin)
%  Output the name of the directory of cache files The basename is in
%  the form of SSSSIIFFFFFF.
%  INPUTS
%   mastcamdata_obj: MASTCAMdata or MASTCAMgroup_wProcCodes that shares same
%   Rover Navigation model, Camera model, and Linearization.
%  OUTPUTS
%   cache_dirname: the name of the directory for the cache files.
%      SSSSIIFFFFFF
%       SSSS: Four-digit sol number after landing day (which was defined as
%             sol 0).
%       II  : Two-digit camera code: “ML” = Mastcam Left (M-34) and 
%             “MR” = Mastcam Right (M-100).
%       FFFFFF: Six-digit sequence number identifier.
%  OPTIONAL Parameters
%   'CAM_CODE' : camera code {'ML', 'MR'}
%     (default) obtained from mastcamdata_obj

if ~isa(mastcamdata_obj,'MASTCAMdata') && ~isa(mastcamdata_obj,'MASTCAMgroup_wProcCode')
    error(['The first input mastcamdata_obj needs to be an instance of ', ...
           'either of MASTCAMdata and MASTCAMgroup_wProcCodes']);
end

cam_code = '';

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'CAM_CODE'
                cam_code = varargin{i+1};
                validatestring(cam_code,{'ML','MR'},mfilename,'CAM_CODE');
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

%%

if iscell(mastcamdata_obj.PRODUCT_ID)
    productID_repre = mastcamdata_obj.PRODUCT_ID{1};
else
    productID_repre = mastcamdata_obj.PRODUCT_ID;
end
propMASTCAMdata = getProp_basenameMASTCAM(productID_repre);

if isempty(cam_code)
    cam_code = propMASTCAMdata.cam_code;
end
if isnumeric(propMASTCAMdata.sol)
    sol = sprintf('%04d',propMASTCAMdata.sol);
end
if isnumeric(propMASTCAMdata.seq_id)
    seq_id = sprintf('%06d',propMASTCAMdata.seq_id);
end
cache_dirname = sprintf('%s%s%s',sol,cam_code,seq_id);

end