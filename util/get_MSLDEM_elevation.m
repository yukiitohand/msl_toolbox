function [elevation_tar] = get_MSLDEM_elevation(northing,easting,MSLDEMdata)

ne_tar = [northing;easting];

L_dem = MSLDEMdata.hdr.lines; S_dem = MSLDEMdata.hdr.samples;
dem_northing   = reshape(MSLDEMdata.hdr.y,L_dem,1);
dem_easting    = reshape(MSLDEMdata.hdr.x,1,S_dem);

yloc = -(northing-dem_northing(1))/MSLDEMdata.hdr.map_info.dy + 1;
xloc = (easting-dem_easting(1))/MSLDEMdata.hdr.map_info.dx + 1;

l_tar = floor(yloc);
s_tar = floor(xloc);

deml_elevation  = MSLDEMdata.lazyEnviReadl(l_tar,0);
demlp1_elevation = MSLDEMdata.lazyEnviReadl(l_tar+1,0);
elevation = cat(1, deml_elevation, demlp1_elevation);
elevation(elevation==MSLDEMdata.hdr.data_ignore_value) = nan;
elevation = elevation(:,[s_tar s_tar+1]);

northeast_l = nan(2,2,2);
northeast_l(1,:,:) = repmat(dem_northing([l_tar l_tar+1])',[1,1,2]);
northeast_l(2,:,:) = repmat(permute(dem_easting([s_tar s_tar+1]),[1,3,2]),[1,2,1]);

for j=1:2
    if j==1
        ppv1 = northeast_l(:,1,1); % plane position vector in image space
        ppv2 = northeast_l(:,1,2);
        ppv3 = northeast_l(:,2,1);
        el1 = elevation(1,1); el2 = elevation(1,2);
        el3 = elevation(2,1);
    elseif j==2
        ppv1 = northeast_l(:,1,2);
        ppv2 = northeast_l(:,2,2);
        ppv3 = northeast_l(:,2,1);
        el1 = elevation(1,2); el2 = elevation(2,2);
        el3 = elevation(2,1);
    end
    [plane_param,is_in_face] = get_plane_param_coefficient(...
            ppv1,ppv2,ppv3,ne_tar,'double',false,false);


    imls_update = find(is_in_face);
    if ~isempty(imls_update)
        elevation_tar(imls_update) = el1 + ([el2 el3] - el1)*plane_param(:,imls_update);
    end
end



end