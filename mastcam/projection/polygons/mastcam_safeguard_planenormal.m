function [mastcam_surfplnc_new,mastcam_emi_new] = mastcam_safeguard_planenormal(...
    mstdata_obj,mastcam_emi,mastcam_surfplnc,rover_nav)

cmmdl_geo = transform_CAHVOR_MODEL_wROVER_NAV(mstdata_obj.CAM_MDL,rover_nav);
cmmdl_geo.get_image_plane_unit_vectors();
cam_C = reshape(cmmdl_geo.C,[1,1,3]);

maxdeg = 90-1;

L_im = mstdata_obj.L_im;
S_im = mstdata_obj.S_im;

[Xc,Yc] = meshgrid(0:(S_im-1),0:(L_im-1));
XYc1d = [Xc(:) Yc(:)]';
pmcc1d = cmmdl_geo.get_p_minus_c_from_xy(XYc1d);
pmcc = reshape(pmcc1d',[L_im,S_im,3]);
lprm = (mastcam_surfplnc(:,:,4) - sum(mastcam_surfplnc(:,:,1:3).*cam_C,3)) ./ sum(pmcc .* mastcam_surfplnc(:,:,1:3),3);
pmcc_vec_intersect = lprm.*pmcc + cam_C;

mastcam_surfplnc_new = mastcam_surfplnc;
mastcam_emi_new      = mastcam_emi     ;

for s=1:S_im
    for l=1:L_im
        if mastcam_emi(l,s)>maxdeg
            % if s==1083 && l==184
            %    keyboard; 
            % end
            n = reshape(mastcam_surfplnc(l,s,1:3),[1,3]);
            pmccls = reshape(pmcc(l,s,:),[1,3]);
            u = cross(n,pmccls);
            u = u ./ sqrt(sum(u.^2));
            theta = deg2rad(mastcam_emi(l,s)-maxdeg);
            Qrot = qGetRotQuaternion( theta , u );
            nrot = qRotatePoint( n, Qrot );
            nrot = reshape(nrot,[1 3]);
            % update plane constant
            pc = nrot * reshape(pmcc_vec_intersect(l,s,:),[3,1]);
            
            % update values 
            mastcam_emi_new(l,s) = maxdeg;
            mastcam_surfplnc_new(l,s,1:3) = reshape(nrot,[1,1,3]);
            mastcam_surfplnc_new(l,s,4) = pc;
            
        end
    end
end


end