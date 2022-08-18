function [mastcam_prj,MSLDEMprj] = proj_mastcam2MSLDEM_v4(mastcamdata_obj,MSLDEMdata,...
    MSLDEMprj,varargin)
% [mastcam_prj] = proj_mastcam2MSLDEM(mastcamdata_obj,MSLDEMdata,...
%     MSLDEMprj,varargin)
%   Project mastcam image onto MSLDEMdata
%  INPUTS:
%    mastcamdata_obj: MASTCAMdata or MASTCAMgroup_eye class object
%    MSLDEMdata: HSI class obj, 
%      MSLDEMdata at
%      https://astrogeology.usgs.gov/search/map/Mars/MarsScienceLaboratory/
%      Mosaics/MSL_Gale_DEM_Mosaic_10m
%   MSLDEMprj: output of the function proj_MSLDEM2mastcam
%    struct having fields:
%    imFOV_mask: boolean image, [L_dem x S_dem x 1]
%      true if in the FOV, false otherwise.
%      FOV should be true if the linear transformation shows the value
%      between -200 and 200+image_edge. The values outside of this are
%      considered to be inaccurate or invalid. In addition, pixels within
%      the distance 50m are considered to be in FOV, just in case.
%    imxy: [[] x []  x 2] 
%      page 1 is the x coordinate in the image frame,
%      page 2 is the y coordinate in the image frame.
%      x and y coordinate values outside of the range between -200 and 
%      200+image_edge are replaced with nans because the values outside of
%      this is far away from the calibration range, therefore accuracy is
%      not guaranteed. Size depend on FOV. The rectangle that minimally
%      encloses the FOV.
%  OUTPUTS:
%    mastcam_prj: struct having fields
%      NorthEastElevation: [L_im x S_im x 3]
%       pages 1,2,3 are northing, easting, and elevation. 
%      MSLDEM_ref: [L_im x S_im x 3]
%       indicate which triangle in the MSLDEMdata, the pixel is located.
%       pages 1,2,3 are sample, line, and j. j is either 1 or 2, indicating
%       which triangle at the (sample, line).
%      range: [L_im x S_im]
%       range for each pixel.
%      MSLDEM_nn: [L_im x S_im x 2]
%       nearest neighbor pixels in DDR image. The first page is the sample
%       indices and the second is the line indexes.
%    MSLDEMprj: one field added to the input MSLDEMprj
%      imFOV_mask_nh: [L_ddr x S_ddr x 1]
%       Boolean, imFOV_mask with hidden points are removed.
%    

is_gpu = false;
precision = 'double';
proc_page = 0;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'PRECISION'
                precision = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

if is_gpu
    % gpu_varargin = {'gpuArray'};
else
    gpu_varargin = {};
end

%-------------------------------------------------------------------------%
% Get the size of the mastcam image
%-------------------------------------------------------------------------%
switch class(mastcamdata_obj)
    case 'MASTCAMdata'
        L_im = mastcamdata_obj.hdr.lines; S_im = mastcamdata_obj.samples;
    case 'MASTCAMgroup_eye'
        L_im = mastcamdata_obj.L_im; S_im = mastcamdata_obj.S_im;
    otherwise
        error('second input needs to be either MASTCAMdata or MASTCAMgroup_eye class.');
end

%-------------------------------------------------------------------------%
% Get cam_C_geo
%-------------------------------------------------------------------------%
cmmdl = mastcamdata_obj.CAM_MDL;
rover_nav_coord = mastcamdata_obj.ROVER_NAV;
cmmdl_geo = transform_CAHVOR_MODEL_wROVER_NAV(cmmdl,rover_nav_coord);
C_geo = reshape(cmmdl_geo.C,3,1);
A_geo = reshape(cmmdl_geo.A,1,3);
% C_geo = reshape(mastcamdata_obj.CAM_MDL_GEO.C_geo,[],1);

switch lower(precision)
    case 'double'
        C_geo = double(C_geo);
    case 'single'
        C_geo = single(C_geo);
end

%-------------------------------------------------------------------------%
% Construct Image grid
%-------------------------------------------------------------------------%
imx_im_1d = 0:(S_im-1);
imy_im_1d = reshape(0:(L_im-1),[],1);
[imx_im,imy_im] = meshgrid(imx_im_1d,imy_im_1d);
im_imxy_vec2d = permute(reshape(cat(2,imx_im,imy_im),[L_im*S_im,2]),[2,1,3]);

PmC = cmmdl_geo.get_p_minus_c_from_xy(im_imxy_vec2d);


%%
%-------------------------------------------------------------------------%
% Get MSL DEM information
%-------------------------------------------------------------------------%
L_dem = MSLDEMdata.hdr.lines; S_dem = MSLDEMdata.hdr.samples;
% dem_northing = reshape(MSLDEMdata.hdr.y(MSLDEMprj.hdr_imxy.line_offset+1):(MSLDEMprj.hdr_imxy.line_offset+MSLDEMprj.hdr_imxy.lines),L_dem,1);
% dem_easting  = reshape(MSLDEMdata.hdr.x,1,S_dem);
% 
% switch lower(precision)
%     case 'double'
%         dem_northing = double(dem_northing);
%         dem_easting = double(dem_easting);
%     case 'single'
%         dem_northing = single(dem_northing);
%         dem_easting = single(dem_easting);
% end

%%
%-------------------------------------------------------------------------%
% Pre-processing before calculation
%-------------------------------------------------------------------------%

% PmC = reshape(mastcamdata_obj.CAM_MDL_GEO.imxy_direc_rov0,[L_im*S_im,3])';

% Placeholders for the image
imxyz_geo = nan(3,L_im*S_im,precision,gpu_varargin{:});
imxyz_geo_ref = nan(3,L_im*S_im,precision,gpu_varargin{:});
imxyz_geo_range = inf(1,L_im*S_im,precision,gpu_varargin{:});
msldem_nn = nan(2,L_im*S_im,precision,gpu_varargin{:});

% if the whole line is outside of FOV, skip
% it is important to take the function any in the dimension 1.
% last line is also removed because of the implementation we did.
% valid_lines = any(MSLDEMprj.imFOV_mask',1);
% valid_lines = find(valid_lines);

% valid_lines = valid_lines(1:end-1);



l1 = MSLDEMprj.hdr_imxy.line_offset+1;
lend = MSLDEMprj.hdr_imxy.line_offset+MSLDEMprj.hdr_imxy.lines;
s1 = MSLDEMprj.hdr_imxy.sample_offset+1;
send = MSLDEMprj.hdr_imxy.sample_offset+MSLDEMprj.hdr_imxy.samples;



MSLDEMdata.fopen_img();
typeName = 'single'; machine = 'ieee-le'; size_type = 4;
len_vs = MSLDEMprj.hdr_imxy.samples;
len_vl = MSLDEMprj.hdr_imxy.lines;
ts_l = size_type*MSLDEMprj.hdr_imxy.sample_offset;
ts_r = size_type*(S_dem-MSLDEMprj.hdr_imxy.sample_offset-len_vs);

fseek(MSLDEMdata.fid_img,size_type*(S_dem*MSLDEMprj.hdr_imxy.line_offset),-1);
fseek(MSLDEMdata.fid_img,ts_l,0);
demlp1_elev = fread(MSLDEMdata.fid_img,len_vs,typeName,0,machine);
fseek(MSLDEMdata.fid_img,ts_r,0);
demlp1_elev(demlp1_elev==MSLDEMdata.hdr.data_ignore_value) = nan;

dem_imFOV_mask_crop = MSLDEMprj.imFOV_mask(l1:lend,s1:send);
dem_imFOVd_mask_crop = ~isnan(MSLDEMprj.imxy(:,:,1));
dem_northing_crop = reshape(MSLDEMdata.hdr.y(l1:lend),[],1);
dem_easting_crop  = reshape(MSLDEMdata.hdr.x(s1:send),1,[]);

% deml_geo = nan(3,len_vs,precision,gpu_varargin{:});
% deml_geo(2,:) = repmat(permute(dem_easting_crop,[1,3,2]),[1,2,1]);

demlp1_geo = nan(3,len_vs,precision,gpu_varargin{:});
demlp1_geo(1,:) = dem_northing_crop(1);
demlp1_geo(2,:) = dem_easting_crop;
demlp1_geo(3,:) = -demlp1_elev;

%figure;

for li = 1:(len_vl-1) % l = 27270 
    %tic;
    l = li+MSLDEMprj.hdr_imxy.line_offset;
    %======================================================================
    % Read elevation of the lth line from MSLDEMdata
    %======================================================================
    % deml_elev  = demlp1_elev;
    fseek(MSLDEMdata.fid_img,ts_l,0);
    demlp1_elev = fread(MSLDEMdata.fid_img,len_vs,typeName,0,machine);
    fseek(MSLDEMdata.fid_img,ts_r,0);
    demlp1_elev(demlp1_elev==MSLDEMdata.hdr.data_ignore_value) = nan;
    
    
    %======================================================================
    % converting to xyz coordinate
    %======================================================================
    % get geographical information from 
    deml_geo = demlp1_geo;
    % deml_geo(1,:) = dem_northing_crop(li);
    % deml_geo(3,:) = -deml_elev;
    
    demlp1_geo(1,:) = dem_northing_crop(li+1);
    demlp1_geo(3,:) = -demlp1_elev;
    
    %======================================================================
    % converting to xyz coordinate
    %====================================================================== 
    deml_imxy = squeeze(MSLDEMprj.imxy(li,:,:))';
    demlp1_imxy = squeeze(MSLDEMprj.imxy(li+1,:,:))';
    
    %======================================================================
    % Evaluate intersecting points of the light of sight and triangles in 
    % the FOV
    %======================================================================
    valid_samples = find(dem_imFOV_mask_crop(li,:));
    if valid_samples(1) > 1
        valid_samples = [valid_samples(1)-1 valid_samples];
    end
    if valid_samples(end) == len_vs
        valid_samples = valid_samples(1:end-1);
    end
    

    for sxy = valid_samples
        for j=1:2
            if j==1
                % plane position vectors in image space
                ppv1 = deml_imxy(:,sxy); 
                ppv2 = deml_imxy(:,sxy+1);
                ppv3 = demlp1_imxy(:,sxy);
                % plane position vectors
                ppv1_geo = deml_geo(:,sxy);
                ppv2_geo = deml_geo(:,sxy+1);
                ppv3_geo = demlp1_geo(:,sxy);
                ppv4_geo = demlp1_geo(:,sxy+1);
                ppv_geo_idxList = [...
                    sxy,sxy+1,  sxy,sxy+1;
                    li,  li,li+1,li+1];
                isinFOV = dem_imFOV_mask_crop(li,sxy).*dem_imFOV_mask_crop(li,sxy+1).*dem_imFOV_mask_crop(li+1,sxy);
                isinFOVd = dem_imFOVd_mask_crop(li,sxy).*dem_imFOVd_mask_crop(li,sxy+1).*dem_imFOVd_mask_crop(li+1,sxy);
            elseif j==2
                % plane position vectors in image space
                ppv1 = deml_imxy(:,sxy+1);
                ppv2 = demlp1_imxy(:,sxy+1);
                ppv3 = demlp1_imxy(:,sxy);
                % plane position vectors
                ppv1_geo = deml_geo(:,sxy+1);
                ppv2_geo = demlp1_geo(:,sxy+1);
                ppv3_geo = demlp1_geo(:,sxy);
                ppv4_geo = deml_geo(:,sxy);
                ppv_geo_idxList = [...
                    sxy+1,sxy+1,  sxy,  sxy;
                      li,li+1,li+1,  li];
                isinFOV = dem_imFOV_mask_crop(li,sxy+1).*dem_imFOV_mask_crop(li+1,sxy+1).*dem_imFOV_mask_crop(li+1,sxy);
                isinFOVd = dem_imFOVd_mask_crop(li,sxy+1).*dem_imFOVd_mask_crop(li+1,sxy+1).*dem_imFOVd_mask_crop(li+1,sxy);
            end
            % ppv_xyList = [ppv1 ppv2 ppv3];
            
            % ppv_geo_idxList = ppv_geo_idxList - ...
            %         [MSLDEMprj.hdr_imxy.sample_offset;...
            %         MSLDEMprj.hdr_imxy.line_offset];
            
            % check if one of the vertices are outside of FOV
            
            if isinFOVd
                xy_min = ceil(min(min(ppv1,ppv2),ppv3));
                xy_max = floor(max(max(ppv1,ppv2),ppv3));
                xy_min = min(max(xy_min,0),[S_im-1;L_im-1]); xy_max = max(min(xy_max,[S_im-1;L_im-1]),0);
                x_list = xy_min(1):xy_max(1); y_list = xy_min(2):xy_max(2);
                idx_xy2d_list = L_im*(x_list) + (y_list+1)';
                idx_xy2d_list = idx_xy2d_list(:);
                [plane_param_im,is_in_face] = get_plane_param_coefficient_2d(...
                    ppv1,ppv2,ppv3,im_imxy_vec2d(:,idx_xy2d_list));
                
            elseif isinFOV
                idx_xy2d_list = 1:(L_im*S_im);
                [line_param,is_intersect] = line_plane_intersect_ldv(...
                        C_geo,PmC,ppv1_geo,ppv2_geo,ppv3_geo,...
                        is_gpu,proc_page);
                is_right_dir = line_param>0;
                pipv = C_geo + PmC(:,is_right_dir).*line_param(:,is_right_dir); % plane intersection position vector
                plane_param = nan(2,L_im*S_im); is_in_face = false(1,L_im*S_im);
                [plane_param(:,is_right_dir),is_in_face(is_right_dir)] = get_plane_param_coefficient_3d(...
                    ppv1_geo,ppv2_geo,ppv3_geo,pipv,precision);
            else
                is_in_face = false;
            end
            if any(is_in_face)
                if isinFOVd
                    [plane_param] = convert_sdtd2st(plane_param_im,...
                         ppv1_geo,ppv2_geo,ppv3_geo,C_geo,A_geo);
                    %plane_param = plane_param_im;
                end
                imxyz_geo_s = ppv1_geo + ([ppv2_geo ppv3_geo] - ppv1_geo)*plane_param;
                imxyz_geo_range_s = sqrt(sum((imxyz_geo_s - C_geo).^2,1));

                imls_update = find(and( is_in_face, imxyz_geo_range_s < imxyz_geo_range(:,idx_xy2d_list) ));
                if ~isempty(imls_update)
                    idxes_update = idx_xy2d_list(imls_update);
                    imxyz_geo(:,idxes_update) = imxyz_geo_s(:,imls_update);
                    imxyz_geo_ref(:,idxes_update) = repmat([sxy;l;j],[1,length(imls_update)]);
                    imxyz_geo_range(idxes_update) = imxyz_geo_range_s(1,imls_update);
                    
                    % get nearest neighbor
                    ppv_geoList = [ppv1_geo ppv2_geo ppv3_geo ppv4_geo];
                    dst = sum((imxyz_geo_s(1:2,imls_update) - permute(ppv_geoList(1:2,:),[1,3,2])).^2,1);
                    [~,min_idx] = min(dst,[],3);
                    msldem_nn(:,idxes_update) = ppv_geo_idxList(:,min_idx);
                end
            end
            %imagesc(reshape(imxyz_geo_range,[L_im,S_im]));
        end
    end
    
    % title(num2str(l));
    % drawnow;
    % toc;
end

MSLDEMdata.fclose_img();

imxyz_geo = permute(reshape(imxyz_geo,[3,L_im,S_im]),[2,3,1]);
imxyz_geo_ref = permute(reshape(imxyz_geo_ref,[3,L_im,S_im]),[2,3,1]);
imxyz_geo_range = reshape(imxyz_geo_range,[L_im,S_im]);
msldem_nn = permute(reshape(msldem_nn,[2,L_im,S_im]),[2,3,1]);

idx_nn = (msldem_nn(:,:,1))*MSLDEMprj.hdr_imxy.lines + msldem_nn(:,:,2);
idx_nn = idx_nn(:);
idx_nn_1d_nisnan = ~isnan(idx_nn);
idx_nn_1d_ok = idx_nn(idx_nn_1d_nisnan);
img_mask_nh = false(MSLDEMprj.hdr_imxy.lines*MSLDEMprj.hdr_imxy.samples,1);
img_mask_nh(idx_nn_1d_ok) = true;
img_mask_nh = reshape(img_mask_nh,[MSLDEMprj.hdr_imxy.lines,MSLDEMprj.hdr_imxy.samples]);

mastcam_prj = [];
mastcam_prj.NorthEastElevation = imxyz_geo;
mastcam_prj.MSLDEM_ref = imxyz_geo_ref;
mastcam_prj.range = imxyz_geo_range;
mastcam_prj.msldem_nn = msldem_nn;

MSLDEMprj.imFOV_mask_nh = img_mask_nh;

end
