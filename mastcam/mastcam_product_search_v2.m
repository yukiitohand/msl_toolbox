function [RDRINDEX_matched] = mastcam_product_search_v2(propMASTCAM)


[volid_ident_support] = mastcam_get_volid_ident_support();
volid_range_mat = cat(1,volid_ident_support.sol_range);

if isnumeric(propMASTCAM.sol)
    sol = propMASTCAM.sol;
    
    volid_candidates = find(and(sol>=volid_range_mat(:,1),sol<=volid_range_mat(:,2)));
    
    basenamePtrn = get_basenameMASTCAM_fromProp(propMASTCAM);
    
    RDRINDEX_matched = [];
    for vi=1:length(volid_candidates)
        volid = volid_candidates(vi);
        load(sprintf('RDRINDEX_MSLMST_%04d_v2.mat',volid),'RDRINDEX_VOLUME_ID',...
        'RDRINDEX_PRODUCT_ID','RDRINDEX_PATH_NAME','RDRINDEX_PLANET_DAY_NUMBER',...
        'RDRINDEX_ROVER_MOTION_COUNTER_SITE','RDRINDEX_ROVER_MOTION_COUNTER_DRIVE',...
        'RDRINDEX_ROVER_MOTION_COUNTER_POSE','RDRINDEX_ROVER_MOTION_COUNTER_RSM',...
        'RDRINDEX_ROVER_MOTION_COUNTER_ARM','RDRINDEX_ROVER_MOTION_COUNTER_CHIMRA',...
        'RDRINDEX_ROVER_MOTION_COUNTER_DRILL');
    
        % first limit the search window to the given sol
        idx_sol = find(RDRINDEX_PLANET_DAY_NUMBER==sol);
        
        matched_result = regexpi(string(RDRINDEX_PRODUCT_ID(idx_sol,:)),basenamePtrn);
        matched_idx = ~isempties(matched_result);
        
        idx_m = idx_sol(matched_idx);
        
        RDRINDEX_matched_vi = struct('VOLUME_ID',cellstr(RDRINDEX_VOLUME_ID(idx_m,:)),...
            'PRODUCT_ID',cellstr(RDRINDEX_PRODUCT_ID(idx_m,:)),...
            'PATH_NAME',cellstr(RDRINDEX_PATH_NAME(idx_m,:)),...
            'ROVER_MOTION_COUNTER_SITE',num2cell(RDRINDEX_ROVER_MOTION_COUNTER_SITE(idx_m)),...
            'ROVER_MOTION_COUNTER_DRIVE',num2cell(RDRINDEX_ROVER_MOTION_COUNTER_DRIVE(idx_m)),...
            'ROVER_MOTION_COUNTER_POSE',num2cell(RDRINDEX_ROVER_MOTION_COUNTER_POSE(idx_m)),...
            'ROVER_MOTION_COUNTER_RSM',num2cell(RDRINDEX_ROVER_MOTION_COUNTER_RSM(idx_m)),...
            'ROVER_MOTION_COUNTER_ARM',num2cell(RDRINDEX_ROVER_MOTION_COUNTER_ARM(idx_m)),...
            'ROVER_MOTION_COUNTER_CHIMRA',num2cell(RDRINDEX_ROVER_MOTION_COUNTER_CHIMRA(idx_m)),...
            'ROVER_MOTION_COUNTER_DRILL',num2cell(RDRINDEX_ROVER_MOTION_COUNTER_DRILL(idx_m)));

        RDRINDEX_matched = [RDRINDEX_matched;RDRINDEX_matched_vi];
       
    end        

end

end