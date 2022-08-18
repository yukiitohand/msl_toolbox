function [mst_rmc] = mastcam_get_RMC(lbl)

[site_id] = mastcam_get_RMC_component(lbl,'SITE');
[drive_id] = mastcam_get_RMC_component(lbl,'DRIVE');
[pose_id] = mastcam_get_RMC_component(lbl,'POSE');
[rsm_mc] = mastcam_get_RMC_component(lbl,'RSM');
[arm_mc] = mastcam_get_RMC_component(lbl,'ARM');
[chimra_mc] = mastcam_get_RMC_component(lbl,'CHIMRA');
[drill_mc] = mastcam_get_RMC_component(lbl,'DRILL');
[hga_mc] = mastcam_get_RMC_component(lbl,'HGA');
[drt_mc] = mastcam_get_RMC_component(lbl,'DRT');
[ic_mc] = mastcam_get_RMC_component(lbl,'IC');


mst_rmc = MSL_RMC('SITE',site_id,'DRIVE',drive_id','POSE',pose_id,...
    'RSM',rsm_mc,'ARM',arm_mc,'CHIMRA',chimra_mc,'DRILL',drill_mc,...
    'HGA',hga_mc,'DRT',drt_mc,'IC',ic_mc);

end

function [rmc_component] = mastcam_get_RMC_component(lbl,component_name)

rmc_component_index = find(strcmpi(component_name,lbl.GROUP_GEOMETRIC_CAMERA_MODEL_PARMS.COORDINATE_SYSTEM_INDEX_NAME));
if isempty(rmc_component_index)
    rmc_component = [];
else
    rmc_component = lbl.GROUP_GEOMETRIC_CAMERA_MODEL_PARMS.REFERENCE_COORD_SYSTEM_INDEX(rmc_component_index);
end

end