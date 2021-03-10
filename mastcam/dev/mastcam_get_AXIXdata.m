function [mstaxixdata] = mastcam_get_AXIXdata(basename,mstgrp_eye,varargin)
% Preample
global msl_env_vars
dirpath_AXIXI = msl_env_vars.dirpath_AXI1;

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            % ## I/O OPTIONS #---------------------------------------------
            case 'DIRPATH'
                dirpath_AXIXI = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

propAXIX = getProp_basenameMASTCAM(basename);
if isempty(regexpi(propAXIX.data_proc_code,'A\dI\d','once'))
    error('Given data is not AXIX data');
end

fpath = joinPath(dirpath_AXIXI,[basename '.IMG']);
if ~exist(fpath,'file')
   fprintf('File %s does not exist.',fpath);
   mstaxixdata = [];
   return;
end

%%
propAXIX_searchtmp = create_propMASTCAMbasename();
propAXIX_search = propAXIX;
propAXIX_search.data_proc_code = propAXIX_searchtmp.data_proc_code;
basename_ref_search = get_basenameMASTCAM_fromProp(propAXIX_search);

active_dpcs = fieldnames(mstgrp_eye.PRODUCT_ID.(propAXIX.product_type)); % active data processing code

mtch_list = [];
for i=1:length(active_dpcs)
    fldnm = active_dpcs{i};
    mtch_list.(fldnm) = find(~isempties(regexpi(mstgrp_eye.PRODUCT_ID.(propAXIX.product_type).(fldnm),...
                           basename_ref_search,'once')));
end
%
if isfield(mtch_list,'DRXX')
    mstdata_ref = mstgrp_eye.(propAXIX.product_type).DRXX(mtch_list.DRXX);
% elseif isfield(mtch_list,'DRXX')
%     mstdata_ref = mstgrp_eye.(propAI.product_type).DRXX(mtch_list.DRXX);
end
mstaxixdata = MASTCAMAIdata(basename,dirpath_AXIXI,'MASTCAMdata_ref',mstdata_ref);


end
