function [mask_hsi,mask_msldemc_UFOVmask_hsiresol,mask_mastcam_hsiresol]...
        = mstproj_get_hsi_mask(MSTproj,hsielem,x_east,y_north,mask_msldemc_UFOVmask)
    mask_hsi = false(hsielem.imszy,hsielem.imszx);
    % nearest HSI image pixel.
    if ~isnan(x_east) && ~isnan(y_north)
        [xi,yi] = hsielem.get_xy_fromNE(x_east,y_north);
        mask_hsi(yi,xi) = true;
    end
    [pos_row,pos_col] = find(mask_msldemc_UFOVmask);
    for i=1:length(pos_row)
        rowi = pos_row(i); coli = pos_col(i);
        [xei,yni] = MSTproj.get_NE_from_msldemc_imUFOVxy(coli,rowi);
        % Get the HSI pixel nearest from one selected pixel in
        % msldemc_UFOVmask.
        [xi,yi] = hsielem.get_xy_fromNE(xei,yni);
        % mask_hsi_ar(i,1) = xi; mask_hsi_ar(i,2) = yi; 
        mask_hsi(yi,xi) = true;
    end
    mask_msldemc_UFOVmask_hsiresol = false(...
        MSTproj.msldemc_imUFOVhdr.lines,...
        MSTproj.msldemc_imUFOVhdr.samples);

    mask_mastcam_hsiresol = false(MSTproj.MASTCAMdata.L_im,...
        MSTproj.MASTCAMdata.S_im);
    [hsirow,hsicol] = find(mask_hsi);

    x1 = MSTproj.msldemc_imUFOVhdr.x(1);
    dx = MSTproj.MSLDEMdata.hdr.map_info.dx;
    y1 = MSTproj.msldemc_imUFOVhdr.y(1);
    dy = MSTproj.MSLDEMdata.hdr.map_info.dy;


    for i=1:length(hsirow)
        hsirowi = hsirow(i); hsicoli = hsicol(i);
        [xrnge_east,yrnge_north] = hsielem.hsi.get_pixel_rangeNE_fromGLTxy(hsicoli,hsirowi);
        xUFOVstrt = ceil((xrnge_east(1)-x1)/dx);
        xUFOVend = floor((xrnge_east(2)-x1)/dx);

        yUFOVstrt = ceil((y1-yrnge_north(2))/dy);
        yUFOVend = floor((y1-yrnge_north(1))/dy);

        for y_demi=yUFOVstrt:yUFOVend
            for x_demi=xUFOVstrt:xUFOVend
                mask_msldemc_UFOVmask_hsiresol(y_demi,x_demi) = true;
                if MSTproj.msldemc_imUFOVmask(y_demi,x_demi)
                    idxes = MSTproj.mapper_msldemc2mastcam_cell{MSTproj.mapper_msldemc2mastcam_mat(y_demi,x_demi)};
                    for ii=1:size(idxes,2)
                        mask_mastcam_hsiresol(idxes(2,ii),idxes(1,ii)) = true;
                        % idxes_jj = obj.MSTproj.mapper_mastcam2msldemc{idxes(2,ii),idxes(1,ii)};
                        % for jj=1:size(idxes_jj,2)
                        %     mask_msldemc_UFOVmask_hsiresol(idxes_jj(2,jj),idxes_jj(1,jj)) = true;
                        % end
                    end
                end
            end
        end

        % [msldemcUFOV_hsii_mask,mask_mastcam_i] = obj.MSTproj.get_rangeUFOVmask(xrnge_east,yrnge_north);
        % mask_msldemc_UFOVmask_hsiresol = or(...
        %     mask_msldemc_UFOVmask_hsiresol,msldemcUFOV_hsii_mask);
        % mask_mastcam_hsiresol = or(...
        %     mask_mastcam_hsiresol,mask_mastcam_i);
    end
end