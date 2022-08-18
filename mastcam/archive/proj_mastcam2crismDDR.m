function [mastcam_prj,DDRprj] = proj_mastcam2crismDDR(mastcamdata_obj,DDRprj,varargin)
% [mastcam_prj,DDRprj] = proj_mastcam2crismDDR(mastcamdata_obj,DDRprj,varargin)
%   Project mastcam image onto crismDDR image
%  INPUTS:
%   mastcamdata_obj: MASTCAMdata or MASTCAMgroup_eye class object
%   DDRprj: struct having fields:
%    imxy: [L_ddr x S_ddr x 2] 
%      page 1 is the x coordinate in the image frame,
%      page 2 is the y coordinate in the image frame.
%      x and y coordinate values outside of the range between -200 and 
%      200+image_edge are replaced with nans because the values outside of
%      this is far away from the calibration range, therefore accuracy is
%      not guaranteed.
%    imFOV_mask: boolean image, [L_ddr x S_ddr x 1]
%      true if in the FOV, false otherwise.
%      FOV should be true if the linear transformation shows the value
%      between -200 and 200+image_edge. The values outside of this are
%      considered to be inaccurate or invalid. In addition, pixels within
%      the distance 50m are considered to be in FOV, just in case.
%    NorthEastElevation: [L_ddr x S_ddr x 3] 
%      pages 1,2,3 are coordinate in north-east-elevation. 
%  OUTPUTS:
%    mastcam_prj: struct having fields
%      NorthEastElevation: [L_im x S_im x 3]
%       pages 1,2,3 are northing, easting, and elevation. 
%      DDR_ref: [L_im x S_im x 3]
%       indicate which triangle in the DDRdata, the pixel is located.
%       pages 1,2,3 are sample, line, and j. j is either 1 or 2, indicating
%       which triangle at the (sample, line).
%      range: [L_im x S_im]
%       range for each pixel.
%      DDR_nn: [L_im x S_im x 2]
%       nearest neighbor pixels in DDR image. The first page is the sample
%       indices and the second is the line indexes.
%    DDRprj: one field added to the input DDRprj
%      imFOV_mask_nh: [L_ddr x S_ddr x 1]
%       Boolean, imFOV_mask with hidden points are removed.

is_gpu = false;
precision = 'double';
proc_page = 0;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'GPU'
                is_gpu = varargin{i+1};
            case 'PRECISION'
                precision = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

% proc_page = batch_size>1;
ddr_geo = DDRprj.NorthEastElevation;
ddr_imFOV_mask = DDRprj.imFOV_mask;
ddr_imxy = DDRprj.imxy;
L_ddr = size(DDRprj.imFOV_mask,1); S_ddr = size(DDRprj.imFOV_mask,2);

cam_C_geo = reshape(mastcamdata_obj.CAM_MDL_GEO.C_geo,3,1);
cam_A_geo = reshape(mastcamdata_obj.CAM_MDL_GEO.A_rov0,1,3);

switch class(mastcamdata_obj)
    case 'MASTCAMdata'
        L_im = mastcamdata_obj.hdr.lines; S_im = mastcamdata_obj.samples;
    case 'MASTCAMgroup_eye'
        L_im = mastcamdata_obj.L_im; S_im = mastcamdata_obj.S_im;
    otherwise
        error('second input needs to be either MASTCAMdata or MASTCAMgroup_eye class.');
end

if is_gpu
    gpu_varargin = {'gpuArray'};
else
    gpu_varargin = {};
end

if is_gpu
    cam_C_geo = gpuArray(cam_C_geo);
    ddr_imFOV_mask = gpuArray(ddr_imFOV_mask);
end

switch lower(precision)
    case 'double'
        cam_C_geo = double(cam_C_geo);
    case 'single'
        cam_C_geo = single(cam_C_geo);
end

% get image grid information
imx_im_1d = 1:S_im;
imy_im_1d = reshape(1:L_im,[],1);
[imx_im,imy_im] = meshgrid(imx_im_1d,imy_im_1d);
imxy_im_2d = permute(reshape(cat(2,imx_im,imy_im),[L_im*S_im,2]),[2,1,3]);

imxy_direc_geo_2d = permute(...
    reshape(mastcamdata_obj.CAM_MDL_GEO.imxy_direc_rov0,[L_im*S_im,3]),...
    [2,1,3]);

% camera center information
cam_C_geo = reshape(cam_C_geo,[],1);

imxyz_geo = nan(3,L_im*S_im,precision,gpu_varargin{:});
imxyz_geo_ref = nan(3,L_im*S_im,precision,gpu_varargin{:});
imxyz_geo_range = inf(1,L_im*S_im,precision,gpu_varargin{:});
ddr_nn = nan(2,L_im*S_im,precision,gpu_varargin{:});

% if the whole line is outside of FOV, skip
% it is important to take the any in the dimension 1.
% last line is also removed because of the implementation we did.
valid_lines = any(ddr_imFOV_mask',1);
valid_lines = find(valid_lines);
valid_lines = valid_lines(1:end-1);
if is_gpu
    valid_lines = gather(valid_lines);
end

ddr_geo(:,:,3) = (-1)*ddr_geo(:,:,3);
ddr_geo = permute(ddr_geo,[3,1,2]);

len_vl = length(valid_lines);
% figure;
for li = 1:len_vl % l = 27270 
    %tic;
    l = valid_lines(li);
    % l
    %==========================================================================
    % projection of ground reference coordinates to ROVER_NAV coordinate
    %==========================================================================
    % get geographical information from 
    ddrl_geo = ddr_geo(:,[l l+1],:);
    
    ddrl_imxy = permute( ddr_imxy([l l+1],:,:),[3,1,2]);

    valid_samples = find(ddr_imFOV_mask(l,:));
    valid_samples = valid_samples(1:end-1);
    % if valid_samples(end) == (hdr_dem_imxy.samples+hdr_dem_imxy.sample_offset)
    %     valid_samples = valid_samples(1:end-1);
    % end

    for s = valid_samples
        
        s_ddr_imxy = s;
        for j=1:2
            if j==1
                ppv1 = ddrl_imxy(:,1,s_ddr_imxy); % plane position vector in image space
                ppv2 = ddrl_imxy(:,1,s_ddr_imxy+1);
                ppv3 = ddrl_imxy(:,2,s_ddr_imxy);
                ppv1_geo = ddrl_geo(:,1,s); % plane position vector
                ppv2_geo = ddrl_geo(:,1,s+1);
                ppv3_geo = ddrl_geo(:,2,s);
                ppv4_geo = ddrl_geo(:,2,s+1);
                ppv_geo_idxList = [...
                    s,s+1,  s,s+1;
                    l,  l,l+1,l+1];
                isinFOV = ddr_imFOV_mask(l,s_ddr_imxy).*ddr_imFOV_mask(l,s_ddr_imxy+1).*ddr_imFOV_mask(l+1,s_ddr_imxy);
            elseif j==2
                ppv1 = ddrl_imxy(:,1,s_ddr_imxy+1);
                ppv2 = ddrl_imxy(:,2,s_ddr_imxy+1);
                ppv3 = ddrl_imxy(:,2,s_ddr_imxy);
                ppv1_geo = ddrl_geo(:,1,s+1);
                ppv2_geo = ddrl_geo(:,2,s+1);
                ppv3_geo = ddrl_geo(:,2,s);
                ppv4_geo = ddrl_geo(:,1,s);
                ppv_geo_idxList = [...
                    s+1,s+1,  s,  s;
                      l,l+1,l+1,  l];
                isinFOV = ddr_imFOV_mask(l,s_ddr_imxy+1).*ddr_imFOV_mask(l+1,s_ddr_imxy+1).*ddr_imFOV_mask(l+1,s_ddr_imxy);
            end
            ppv_xyList = [ppv1 ppv2 ppv3];
            ppv_geoList = [ppv1_geo ppv2_geo ppv3_geo ppv4_geo];
            
            if any(isnan(ppv_xyList(:)))
                ppv_xy_isnan = true;
            else
                ppv_xy_isnan = false;
            end 
            % test if the intersection is within the triangle or not.
            % first evaluate which one to use
            %if ppv_xy_isnan, cnd_xy = nan;
            %else
            %    cnd_xy = cond([ppv2 ppv3] - ppv1);
            %end
            %cnd_geo = cond([ppv2_geo ppv3_geo] - ppv1_geo);
            if ~ppv_xy_isnan % && cnd_xy < cnd_geo*10 && 
                % if (l==414) && (s==209) && (j==2)
                %    1;
                % end
                xy_min = ceil(min(ppv_xyList,[],2));
                xy_max = floor(max(ppv_xyList,[],2));

                xy_min = min(max(xy_min,[1;1]),[S_im;L_im]); xy_max = max(min(xy_max,[S_im;L_im]),[1;1]);
                x_list = xy_min(1):xy_max(1); y_list = xy_min(2):xy_max(2);
                idx_xy2d_list = L_im*(x_list-1) + y_list';
                idx_xy2d_list = idx_xy2d_list(:);
                [plane_param,is_in_face] = get_plane_param_coefficient(...
                    ppv1,ppv2,ppv3,imxy_im_2d(:,idx_xy2d_list),precision,is_gpu,proc_page);
                % [plane_param] = convert_sdtd2st(plane_param_im,...
                %    ppv1_geo,ppv2_geo,ppv3_geo,cam_C_geo,cam_A_geo);
            else
                if ~isinFOV
                    is_in_face = false;
                else
                    % test line segment intersect with the plane determined by the
                    % three points
                    % if (l==414) && (s==209) && (j==2)
                    %     1;
                    % end
                    idx_xy2d_list = 1:(L_im*S_im);
                    [line_param,is_intersect] = line_plane_intersect_ldv(...
                            cam_C_geo,imxy_direc_geo_2d,ppv1_geo,ppv2_geo,ppv3_geo,...
                            is_gpu,proc_page);
                    is_right_dir = line_param>0;
                    pipv = cam_C_geo + imxy_direc_geo_2d(:,is_right_dir).*line_param(:,is_right_dir); % plane intersection position vector
                    plane_param = nan(2,L_im*S_im);
                    is_in_face = false(1,L_im*S_im);
                    [plane_param(:,is_right_dir),is_in_face(is_right_dir)] = get_plane_param_coefficient(...
                        ppv1_geo,ppv2_geo,ppv3_geo,pipv,precision,is_gpu,proc_page);
                end
            end
            
            if any(is_in_face)
                imxyz_geo_s = ppv1_geo + ([ppv2_geo ppv3_geo] - ppv1_geo) *plane_param;
                imxyz_geo_range_s = sqrt(sum((imxyz_geo_s - cam_C_geo).^2,1));

                imls_update = find(and( is_in_face, imxyz_geo_range_s < imxyz_geo_range(:,idx_xy2d_list) ));
                if ~isempty(imls_update)
                    idxes_update = idx_xy2d_list(imls_update);
                    imxyz_geo(:,idxes_update) = imxyz_geo_s(:,imls_update);
                    imxyz_geo_ref(:,idxes_update) = repmat([s;l;j],[1,length(imls_update)]);
                    imxyz_geo_range(idxes_update) = imxyz_geo_range_s(1,imls_update);
                    
                    % get nearest neighbor
                    dst = sum((imxyz_geo_s(1:2,imls_update) - permute(ppv_geoList(1:2,:),[1,3,2])).^2,1);
                    [~,min_idx] = min(dst,[],3);
                    ddr_nn(:,idxes_update) = ppv_geo_idxList(:,min_idx);
                end
            end
        end
    end
    %toc;
    % imagesc(reshape(imxyz_geo_range,[L_im,S_im]));
    % title(num2str(l));
    % drawnow;
end

if is_gpu
    [imxyz_geo,imxyz_geo_ref,imxyz_geo_range] = gather(imxyz_geo,imxyz_geo_ref,imxyz_geo_range);

end

ddr_nn = permute(reshape(ddr_nn,[2,L_im,S_im]),[2,3,1]);
% compute which pixels are in the MASTCAM image. In particualr, this
% evaluates whether or not a pixel is in the image in a nearest neighbor
% sense. 
% idx_nn = reshape(ddr_nn,[L_im*S_im,2]);
idx_nn = (ddr_nn(:,:,1)-1)*L_ddr + ddr_nn(:,:,2);
idx_nn = idx_nn(:);
idx_nn_1d_nisnan = ~isnan(idx_nn);
idx_nn_1d_ok = idx_nn(idx_nn_1d_nisnan);
img_mask_nh = false(L_ddr*S_ddr,1);
img_mask_nh(idx_nn_1d_ok) = true;
img_mask_nh = reshape(img_mask_nh,L_ddr,S_ddr);


mastcam_prj = [];
mastcam_prj.NorthEastElevation = permute(reshape(imxyz_geo,[3,L_im,S_im]),[2,3,1]);
mastcam_prj.DDR_ref = permute(reshape(imxyz_geo_ref,[3,L_im,S_im]),[2,3,1]);
mastcam_prj.range = reshape(imxyz_geo_range,[L_im,S_im]);
mastcam_prj.DDR_nn = ddr_nn;

DDRprj.imFOV_mask_nh = img_mask_nh;

end
