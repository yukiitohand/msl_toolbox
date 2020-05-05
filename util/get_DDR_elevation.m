function [elevation_tar,ddr_ref,ddr_nn] = get_DDR_elevation(northing,easting,crism_DDRdata)

[nee] = get_CRISM_NorthEastElevation(crism_DDRdata,'RE',3396190);

northeast = nee(:,:,[1,2]);
elevation = nee(:,:,3);

ne_tar = [northing;easting];
elevation_tar = nan(1,size(ne_tar,2));
ddr_ref = nan(3,size(ne_tar,2));
ddr_nn = nan(2,size(ne_tar,2));

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
                ppv4 = northeast_l(:,2,s+1);
                el1 = elevation(l,s); el2 = elevation(l,s+1);
                el3 = elevation(l+1,s);
                ppv_idxList = [...
                    l,  l,l+1,l+1;
                    s,s+1,  s,s+1];
            elseif j==2
                ppv1 = northeast_l(:,1,s+1);
                ppv2 = northeast_l(:,2,s+1);
                ppv3 = northeast_l(:,2,s);
                ppv4 = northeast_l(:,1,s);
                el1 = elevation(l,s+1); el2 = elevation(l+1,s+1);
                el3 = elevation(l+1,s);
                ppv_idxList = [...
                      l,l+1,l+1,  l;
                    s+1,s+1,  s,  s];
            end
            ppv_List = [ppv1 ppv2 ppv3 ppv4];
            
            [plane_param,is_in_face] = get_plane_param_coefficient(...
                    ppv1,ppv2,ppv3,ne_tar,'double',false,false);
            
            imls_update = find(is_in_face);
            if ~isempty(imls_update)
                imxy_s = ppv1 + ([ppv2 ppv3] - ppv1) *plane_param(:,imls_update);
                elevation_tar(imls_update) = el1 + ([el2 el3] - el1)*plane_param(:,imls_update);
                
                ddr_ref(:,imls_update) = repmat([s;l;j],[1,length(imls_update)]);
                % get nearest neighbor
                dst = sum((imxy_s(1:2,imls_update) - permute(ppv_List(1:2,:),[1,3,2])).^2,1);
                [~,min_idx] = min(dst,[],3);
                ddr_nn(:,imls_update) = ppv_idxList(:,min_idx);
            end
        end
    end
    toc;
end


end
