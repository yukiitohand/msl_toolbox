function [RDRINDEX_matched] = mastcam_product_search(propMASTCAM)


[volid_ident_support] = mastcam_get_volid_ident_support();
volid_range_mat = cat(1,volid_ident_support.sol_range);

if isnumeric(propMASTCAM.sol)
    sol = propMASTCAM.sol;
    
    volid_candidates = find(and(sol>=volid_range_mat(:,1),sol<=volid_range_mat(:,2)));
    
    basenamePtrn = get_basenameMASTCAM_fromProp(propMASTCAM);
    
    RDRINDEX_matched = [];
    for vi=1:length(volid_candidates)
        volid = volid_candidates(vi);
        load(sprintf('RDRINDEX_MSLMST_%04d.mat',volid),'rdrindex');
        [m,mi] = searchby('PRODUCT_ID',basenamePtrn,rdrindex);
        RDRINDEX_matched = [RDRINDEX_matched; m];
    end
    
    
    
end