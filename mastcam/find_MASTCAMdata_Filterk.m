function [mstdatak] = find_MASTCAMdata_Filterk(mstgrp_wpc,k,priority_dtype_list)

mstdatak = [];
j=1; Ldl = length(priority_dtype_list);
while isempty(mstdatak) &&  ( j <= Ldl )
    dtype = priority_dtype_list{j};
    if ~isempty(mstgrp_wpc.(dtype))
        ii=1; Ld = length(mstgrp_wpc.(dtype));
        while isempty(mstdatak) && ( ii <= Ld )
            if mstgrp_wpc.(dtype)(ii).FILTER_NUMBER==k
                mstdatak = mstgrp_wpc.(dtype)(ii);
            end
            ii = ii+1;
        end
    end
    j = j+1;
end

end