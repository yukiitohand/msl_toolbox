function [msldem_imFOVmask] = mastcam_get_imFOVmask_msldem_cahvor_lr1(...
    MSLDEMdata,L_im,S_im,cmmdl_geo,coef_mrgn,proc_mode)

cmmdl_geo.get_image_plane_unit_vectors();

% 
% fprintf('%f,%f,%f\n',cmmdl_geo.C(1),cmmdl_geo.C(2),cmmdl_geo.C(3));
switch upper(proc_mode)
    case 'ENCLOSING_RECTANGLE'
        tic; [msldem_imFOVmask] = cahv_get_imFOVmask_MSLDEM_enclosing_rectangular_mex(...
            MSLDEMdata.imgpath,MSLDEMdata.hdr,MSLDEMdata.hdr.y,MSLDEMdata.hdr.x,...
            S_im,L_im,cmmdl_geo,coef_mrgn); toc;
    case 'SURROUNDING_COMPLEMENT'
        tic; 
        [msldem_imFOVmask] = cahv_get_imFOVmask_MSLDEM_surrounding_complement_fast_mex(...
            MSLDEMdata.imgpath,MSLDEMdata.hdr,MSLDEMdata.hdr.y,MSLDEMdata.hdr.x,...
            S_im,L_im,cmmdl_geo,coef_mrgn);
        % msldem_imFOVmask(msldem_imFOVmask<0) = 0;
        toc;
    otherwise
        error('Undefined imFOVmask_mode %s.',proc_mode);
end




end