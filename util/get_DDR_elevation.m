function [elevation_tar] = get_DDR_elevation(northing,easting,DDRprj)

northeast = DDRprj.NorthEastElevation(:,:,[1,2]);

elevation = DDRprj.NorthEastElevation(:,:,3);

ne_tar = [northing;easting];
elevation_tar = nan(1,size(ne_tar,2));

L_ddr = size(northeast,1); S_ddr = size(northeast,2);
for l = 1:(L_ddr-1)
    tic;
    %==========================================================================
    % projection of ground reference coordinates to ROVER_NAV coordinate
    %==========================================================================
    
    northeast_l = permute( northeast([l l+1],:,:),[3,1,2]);
    for s = 1:(S_ddr-1)
        for j=1:2
            if j==1
                ppv1 = northeast_l(:,1,s); % plane position vector in image space
                ppv2 = northeast_l(:,1,s+1);
                ppv3 = northeast_l(:,2,s);
                el1 = elevation(l,s); el2 = elevation(l,s+1);
                el3 = elevation(l+1,s);
            elseif j==2
                ppv1 = northeast_l(:,1,s+1);
                ppv2 = northeast_l(:,2,s+1);
                ppv3 = northeast_l(:,2,s);
                el1 = elevation(l,s+1); el2 = elevation(l+1,s+1);
                el3 = elevation(l+1,s);
            end
            [plane_param,is_in_face] = get_plane_param_coefficient(...
                    ppv1,ppv2,ppv3,ne_tar,'double',false,false);
            
            
            imls_update = find(is_in_face);
            if ~isempty(imls_update)
                elevation_tar(imls_update) = el1 + ([el2 el3] - el1)*plane_param(:,imls_update);
            end
        end
    end
    toc;
end


end
