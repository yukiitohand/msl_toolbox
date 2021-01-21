function [msldem_imFOVmask,lrange,srange] = mastcam_get_imFOVmask_msldem_cahvor_lr0(...
    MSLDEMdata,L_im,S_im,cmmdl_geo,coef_mrgn,proc_mode)

cmmdl_geo.get_image_plane_unit_vectors();

% define permeter pixels.
perimeter_left = cat(1,-0.5*ones(1,L_im+1),(-0.5):(L_im-0.5));
perimeter_right = cat(1,(S_im-0.5)*ones(1,L_im+1),(-0.5):(L_im-0.5));
perimeter_upper = cat(1,(-0.5):(S_im-0.5),-0.5*ones(1,S_im+1));
perimeter_btm = cat(1,(-0.5):(S_im-0.5),(L_im-0.5)*ones(1,S_im+1));

% get apparent permieter pixel coordinate values
xy_ap_pl = cmmdl_geo.get_apparent_xy_from_xy(perimeter_left);
xy_ap_pr = cmmdl_geo.get_apparent_xy_from_xy(perimeter_right);
xy_ap_pu = cmmdl_geo.get_apparent_xy_from_xy(perimeter_upper);
xy_ap_pb = cmmdl_geo.get_apparent_xy_from_xy(perimeter_btm);

% get the range
srange = [min(xy_ap_pl(1,:),[],2),max(xy_ap_pr(1,:),[],2)];
lrange = [min(xy_ap_pu(2,:),[],2),max(xy_ap_pb(2,:),[],2)];

% 
% fprintf('%f,%f,%f\n',cmmdl_geo.C(1),cmmdl_geo.C(2),cmmdl_geo.C(3));
switch upper(proc_mode)
    case 'ENCLOSING_RECTANGLE'
        tic; [msldem_imFOVmask] = cahvor_get_imFOVmask_MSLDEM_enclosing_rectangle_mex(...
            MSLDEMdata.imgpath,MSLDEMdata.hdr,MSLDEMdata.northing,MSLDEMdata.easting,...
            srange,lrange,cmmdl_geo,coef_mrgn); toc;

    case 'SURROUNDING_COMPLEMENT'
        tic;
        [msldem_imFOVmask] = cahvor_get_imFOVmask_MSLDEM_surrounding_complement_fast_mex(...
            MSLDEMdata.imgpath,MSLDEMdata.hdr,MSLDEMdata.northing,MSLDEMdata.easting,...
            S_im,L_im,srange,lrange,cmmdl_geo,coef_mrgn);
        % msldem_imFOVmask(msldem_imFOVmask<0) = 0;
        toc;
        
    otherwise
        error('Undefined imFOVmask_mode %s.',proc_mode);
end




end
        
        



