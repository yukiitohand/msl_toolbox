function [msldem_imFOVmask] = mastcam_get_imFOVmask_msldem_cahvor_iaumars_lr1(...
    MSLDEMdata,L_im,S_im,cmmdl_iaumars,coef_mrgn,proc_mode,mars_re,mars_rp)

cmmdl_iaumars.get_image_plane_unit_vectors();

mslrad_latitude = deg2rad(MSLDEMdata.latitude);
mslrad_longitude = deg2rad(MSLDEMdata.longitude);

% 
% fprintf('%f,%f,%f\n',cmmdl_geo.C(1),cmmdl_geo.C(2),cmmdl_geo.C(3));
switch upper(proc_mode)
    case 'ENCLOSING_RECTANGLE'
        error('Not implemented yet');
    case 'SURROUNDING_COMPLEMENT'
        tic;
        [msldem_imFOVmask] = cahv_iaumars_get_imFOVmask_MSLDEM_scf2_mex(...
            MSLDEMdata.imgpath,MSLDEMdata.hdr,mslrad_latitude,mslrad_longitude,...
            MSLDEMdata.OFFSET,S_im,L_im,cmmdl_iaumars,coef_mrgn,mars_re,mars_rp);
        % msldem_imFOVmask(msldem_imFOVmask<0) = 0;
        toc;
    otherwise
        error('Undefined imFOVmask_mode %s.',proc_mode);
end




end