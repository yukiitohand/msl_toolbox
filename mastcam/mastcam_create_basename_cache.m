function [basename_cache_com] = mastcam_create_basename_cache(mastcamdata_obj)

if iscell(mastcamdata_obj.PRODUCT_ID)
    productID_repre = mastcamdata_obj.PRODUCT_ID{1};
else
    productID_repre = mastcamdata_obj.PRODUCT_ID;
end
propMASTCAMdata = getProp_basenameMASTCAM(productID_repre);
cam_code = propMASTCAMdata.cam_code;
if isnumeric(propMASTCAMdata.sol)
    sol = sprintf('%04d',propMASTCAMdata.sol);
end
if isnumeric(propMASTCAMdata.seq_id)
    seq_id = sprintf('%06d',propMASTCAMdata.seq_id);
end
site_id  = mastcamdata_obj.RMC.SITE;
drive_id = mastcamdata_obj.RMC.DRIVE;
pose_id  = mastcamdata_obj.RMC.POSE;
rsm_mc   = mastcamdata_obj.RMC.RSM;
basename_cache_com = sprintf('%s%s%s_SITE%03dDRIVE%04dPOSE%03dRSM%03d',sol,cam_code,seq_id,...
    site_id,drive_id,pose_id,rsm_mc);

end