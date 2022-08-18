function [ifovH,ifovV] = mastcam_get_pixel_resolution(mstdata_obj,...
    mastcam_range,mastcam_emi,mastcam_surfplnc,rover_nav)

cmmdl_geo = transform_CAHVOR_MODEL_wROVER_NAV(mstdata_obj.CAM_MDL,rover_nav);
cmmdl_geo.get_image_plane_unit_vectors();

L_im = mstdata_obj.L_im;
S_im = mstdata_obj.S_im;

% hd = cmmdl_geo.Hdash;
% vd = cmmdl_geo.Vdash;
% hd_horiz = sqrt(sum(hd(1:2).^2));
% hd_vert = hd(3);
% 
% vd_horiz = sqrt(sum(vd(1:2).^2));
% vd_vert = vd(3);

cam_C = reshape(cmmdl_geo.C,[1,1,3]);

[Xl,Yl] = meshgrid(0:(S_im-1),-0.5:(L_im-0.5));
[Xs,Ys] = meshgrid(-0.5:(S_im-0.5),0:(L_im-1));
% [Xc,Yc] = meshgrid(0:(S_im-1),0:(L_im-1));

XYl1d = [Xl(:) Yl(:)]';
pmcl1d = cmmdl_geo.get_p_minus_c_from_xy(XYl1d);
pmcl = reshape(pmcl1d',[L_im+1,S_im,3]);
% pmcl = pmcl ./ sqrt(sum(pmcl.^2,3));
% ifovvert = abs(acos(sum(pmcl(1:end-1,:,:).*pmcl(2:end,:,:),3)));

% Evaluate the intersection of the pmcl with each surface plane
lprmu = (mastcam_surfplnc(:,:,4) - sum(mastcam_surfplnc(:,:,1:3).*cam_C,3)) ./ sum(pmcl(1:end-1,:,:) .* mastcam_surfplnc(:,:,1:3),3);
lprmd = (mastcam_surfplnc(:,:,4) - sum(mastcam_surfplnc(:,:,1:3).*cam_C,3)) ./ sum(pmcl(2:end,:,:) .* mastcam_surfplnc(:,:,1:3),3);

ifovV = sqrt(sum((lprmu.*pmcl(1:end-1,:,:) - lprmd.*pmcl(2:end,:,:)).^2,3));


XYs1d = [Xs(:) Ys(:)]';
pmcs1d = cmmdl_geo.get_p_minus_c_from_xy(XYs1d);
pmcs = reshape(pmcs1d',[L_im,S_im+1,3]);
% pmcs = pmcs ./ sqrt(sum(pmcs.^2,3));
% ifovhoriz = abs(acos(sum(pmcs(:,1:end-1,:).*pmcs(:,2:end,:),3)));

lprml = (mastcam_surfplnc(:,:,4) - sum(mastcam_surfplnc(:,:,1:3).*cam_C,3)) ./ sum(pmcs(:,1:end-1,:) .* mastcam_surfplnc(:,:,1:3),3);
lprmr = (mastcam_surfplnc(:,:,4) - sum(mastcam_surfplnc(:,:,1:3).*cam_C,3)) ./ sum(pmcs(:,2:end,:) .* mastcam_surfplnc(:,:,1:3),3);
ifovH = sqrt(sum((lprml.*pmcs(:,1:end-1,:) - lprmr.*pmcs(:,2:end,:)).^2,3));


% XYc1d = [Xc(:) Yc(:)]';
% pmcc1d = cmmdl_geo.get_p_minus_c_from_xy(XYc1d);
% pmcc = reshape(pmcc1d',[L_im,S_im,3]);
% pmcc = pmcc ./ sqrt(sum(pmcc.^2,3));
% % pmcc_ang = cat(3,sqrt(sum(pmcc(:,:,1:2).^2,2)),abs(pmcc(:,:,3)));
% % theta = atan(abs(pmcc(:,:,3))./sqrt(sum(pmcc(:,:,1:2).^2,3)));
% 
% ifovH = 2 * mastcam_range .* ifovhoriz .* sqrt(hd_horiz.^2 + (hd_vert ./ sin(theta)).^2);
% ifovV = 2 * mastcam_range .* ifovvert  .* sqrt(vd_horiz.^2 + (vd_vert ./ sin(theta)).^2);

end
