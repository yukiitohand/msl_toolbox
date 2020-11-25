function [mstaxixdata] = mastcam_get_AXIXdata_from_ref(mstdata_ref,varargin)
% Preample
global msl_env_vars
dirpath_AXIX = msl_env_vars.dirpath_AXI1;
vr_AXIX = 1;

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            % ## I/O OPTIONS #---------------------------------------------
            case 'DIRPATH'
                dirpath_AXIX = varargin{i+1};
            case 'VEERSION'
                vr_AXIX = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

propAXIX_search = getProp_basenameMASTCAM(mstdata_ref.basename);
propAXIX_search.data_proc_code = sprintf('A%1dI%1d',mstdata_ref.FILTER_NUMBER,vr_AXIX);
basenameAXIX_search = get_basenameMASTCAM_fromProp(propAXIX_search);

list_dir = dir(dirpath_AXIX);
fnamelist = {list_dir.name};
basenameAXIX = extractMatchedBasename_v2(basenameAXIX_search,fnamelist);

fpath_AXIX = joinPath(dirpath_AXIX,[basenameAXIX '.IMG']);
if ~exist(fpath_AXIX,'file')
   fprintf('File %s does not exist.',fpath_AXIX);
   mstaxixdata = [];
   return;
end

mstaxixdata = MASTCAMdataAXIX(basenameAXIX,dirpath_AXIX,'MASTCAMdata_ref',mstdata_ref);


end