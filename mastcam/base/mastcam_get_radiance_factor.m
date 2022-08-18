function [rf,rofst] = mastcam_get_radiance_factor(mst_instID,filter_number,lbl)

rf_raw = lbl.GROUP_PROCESSING_PARMS.RADIANCE_SCALING_FACTOR;
rofst_raw = lbl.GROUP_PROCESSING_PARMS.RADIANCE_OFFSET;

switch upper(mst_instID)
    case {'MAST_LEFT','L'}
        switch filter_number
            case 0
                rf = rf_raw;    rofst  = rofst_raw;
            case 1
                rf = rf_raw(2); rofst = rofst_raw(2);
            case 2
                rf = rf_raw(3); rofst = rofst_raw(3);
            case 3
                rf = rf_raw(1); rofst = rofst_raw(1);
            case 4
                rf = rf_raw(1); rofst = rofst_raw(1);
            case 5
                % should be same for all elements
                rf = rf_raw(1); rofst = rofst_raw(1);
            case 6
                rf = rf_raw(1); rofst = rofst_raw(1);
        end

    case {'MAST_RIGHT','R'}
        switch filter_number
            case 0
                rf = rf_raw;    rofst = rofst_raw;
            case 1
                rf = rf_raw(2); rofst = rofst_raw(2);
            case 2
                rf = rf_raw(3); rofst = rofst_raw(3);
            case 3
                rf = rf_raw(1); rofst = rofst_raw(1);
            case 4
                % should be same for all elements
                rf = rf_raw(1); rofst = rofst_raw(1);
            case 5
                % should be same for all elements
                rf = rf_raw(1); rofst = rofst_raw(1);
            case 6
                % should be same for all elements
                rf = rf_raw(1); rofst = rofst_raw(1);
        end
    otherwise
        error('Undefined Eye %s. Either {"L","R"} is accepted.',mst_instID);
end