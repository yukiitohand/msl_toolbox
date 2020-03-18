function [mst_cahvor_model] = mastcam_get_cahvor_model(lbl)

[cmmdl_C] = mastcam_get_cahvor_model_component(lbl,'C');
[cmmdl_A] = mastcam_get_cahvor_model_component(lbl,'A');
[cmmdl_H] = mastcam_get_cahvor_model_component(lbl,'H');
[cmmdl_V] = mastcam_get_cahvor_model_component(lbl,'V');
[cmmdl_O] = mastcam_get_cahvor_model_component(lbl,'O');
[cmmdl_R] = mastcam_get_cahvor_model_component(lbl,'R');

mst_cahvor_model = CAHVOR_MODEL('C',cmmdl_C,'A',cmmdl_A,'H',cmmdl_H,...
    'V',cmmdl_V,'O',cmmdl_O,'R',cmmdl_R);

end

function [cmmdl_component] = mastcam_get_cahvor_model_component(lbl,component_id)

cmmdl_component_index = find(strcmpi(component_id, lbl.GROUP_GEOMETRIC_CAMERA_MODEL_PARMS.MODEL_COMPONENT_ID));
if isempty(cmmdl_component_index)
    cmmdl_component = [];
else
    cmmdl_component = lbl.GROUP_GEOMETRIC_CAMERA_MODEL_PARMS.(sprintf('MODEL_COMPONENT_%1d',cmmdl_component_index));
end

end