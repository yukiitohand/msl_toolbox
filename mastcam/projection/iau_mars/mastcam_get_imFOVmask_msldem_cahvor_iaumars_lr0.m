function [msldem_imFOVmask,lrange,srange] = mastcam_get_imFOVmask_msldem_cahvor_iaumars_lr0(...
    MSLDEMdata,L_im,S_im,cmmdl_iaumars,coef_mrgn,proc_mode,mars_re,mars_rp)
% Wrapper function for the 

cmmdl_iaumars.get_image_plane_unit_vectors();

% define permeter pixels.
perimeter_left = cat(1,-0.5*ones(1,L_im+1),(-0.5):(L_im-0.5));
perimeter_right = cat(1,(S_im-0.5)*ones(1,L_im+1),(-0.5):(L_im-0.5));
perimeter_upper = cat(1,(-0.5):(S_im-0.5),-0.5*ones(1,S_im+1));
perimeter_btm = cat(1,(-0.5):(S_im-0.5),(L_im-0.5)*ones(1,S_im+1));

% get apparent permieter pixel coordinate values
xy_ap_pl = cmmdl_iaumars.get_apparent_xy_from_xy(perimeter_left);
xy_ap_pr = cmmdl_iaumars.get_apparent_xy_from_xy(perimeter_right);
xy_ap_pu = cmmdl_iaumars.get_apparent_xy_from_xy(perimeter_upper);
xy_ap_pb = cmmdl_iaumars.get_apparent_xy_from_xy(perimeter_btm);

% get the range
srange = [min(xy_ap_pl(1,:),[],2),max(xy_ap_pr(1,:),[],2)];
lrange = [min(xy_ap_pu(2,:),[],2),max(xy_ap_pb(2,:),[],2)];

mslrad_latitude = deg2rad(MSLDEMdata.latitude);
mslrad_longitude = deg2rad(MSLDEMdata.longitude);

% 
% fprintf('%f,%f,%f\n',cmmdl_geo.C(1),cmmdl_geo.C(2),cmmdl_geo.C(3));
switch upper(proc_mode)
    case 'ENCLOSING_RECTANGLE'
       error('Not implemented yet');

    case 'SURROUNDING_COMPLEMENT'
        tic;
        [msldem_imFOVmask] = cahvor_iaumars_get_imFOVmask_MSLDEM_scf2_mex(...
            MSLDEMdata.imgpath,MSLDEMdata.hdr, ...
            mslrad_latitude,mslrad_longitude,MSLDEMdata.OFFSET, ...
            S_im,L_im,srange,lrange,cmmdl_iaumars,coef_mrgn,mars_re,mars_rp);
        % msldem_imFOVmask(msldem_imFOVmask<0) = 0;
        toc;
        
    otherwise
        error('Undefined imFOVmask_mode %s.',proc_mode);
end




end
     