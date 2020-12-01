classdef MASTCAMCameraProjectionMSLDEM < handle
    %MASTCAMCameraProjectionMSLDEM
    %  class handling camera projection.
    %  Properties
    %    MASTCAMdata: MASTCAMdata or MASTCAMgroup_eye class object
    %    MSLDEMdata: HSI class obj, 
    %      MSLDEMdata at
    %      https://astrogeology.usgs.gov/search/map/Mars/MarsScienceLaboratory/
    %      Mosaics/MSL_Gale_DEM_Mosaic_10m
    %    msldemc_imFOVmask: boolean image, [L_demc x S_demc x 1]
    %      true if in the FOV, false otherwise.
    %      FOV should be true if the linear transformation shows the value
    %      between -200 and 200+image_edge. The values outside of this are
    %      considered to be inaccurate or invalid. In addition, pixels within
    %      the distance 50m are considered to be in FOV, just in case.
    %    msldemc_imFOVxy: [L_demc x L_demc x 2] 
    %      page 1 is the x coordinate in the image frame,
    %      page 2 is the y coordinate in the image frame.
    %      x and y coordinate values outside of the range between -200 and 
    %      200+image_edge are replaced with nans because the values outside of
    %      this is far away from the calibration range, therefore accuracy is
    %      not guaranteed. Size depend on FOV. The rectangle that minimally
    %      encloses the FOV.
    %    mastcam_NEE: [L_im x S_im x 3]
    %      pages 1,2,3 are northing, easting, and elevation. 
    %    mastcam_msldemc_ref: [L_im x S_im x 3]
    %      indicate which triangle in the MSLDEMdata, the pixel is located.
    %      pages 1,2,3 are sample, line, and j. j is either 1 or 2, indicating
    %      which triangle at the (sample, line).
    %     mastcam_range: [L_im x S_im]
    %      range for each pixel.
    %     mastcam_msldemc_nn: [L_im x S_im x 2]
    %      nearest neighbor pixels in DDR image. The first page is the sample
    %      indices and the second is the line indexes.
    %     msldemc_imFOVmask_ctrnn: [L_demc x S_demc x 1]
    %       Boolean, imFOV_mask with hidden points are removed.
    %     mastcam_emi: [L_im x S_im x 1]
    %       emission angle in degree
    %     mastcam_surfplnc: [L_im x S_im x 4]
    %       parameters for the surface plane equation. Surface plane is
    %       defined as the triangle that intersects with central pixel pointing
    %       vectors. The first three layers corresponds to 
    %     
    
    properties
        MASTCAMdata
        MSLDEMdata
        FOVcone
        msldemc_imFOVhdr
        msldemc_imFOVmask
        msldemc_imFOVxy
        msldemc_imFOVmask_ctrnn % 1 for nearest DEM pixels for the center of image pixels.
        mastcam_NEE
        mastcam_range
        mastcam_msldemc_ref
        mastcam_msldemc_nn
        qNEEpolygons          % quick NEE polygons
        qmastcam_imxypolygons % quick imxy polygons
        msldemc_imUFOVhdr  % offset information of unobstructed FOV
        msldemc_imUFOVmask % unobstructed FOV mask
        msldemc_imUFOVxy   % camera image (x,y) in the UFOV 
        msldemc_imUFOVxynn % nearest integer (x,y) pixels
        mastcam_msldemc_nn_UFOVmask
        mastcam_emi
        mastcam_surfplnc
        mapper_mastcam2msldemc
        mapper_msldemc2mastcam_cell
        mapper_msldemc2mastcam_mat
        cachedirpath
        cachefilepath
        cachebasename_common
        version
        rover_nav_version
    end
    
    methods
        function obj = MASTCAMCameraProjectionMSLDEM(mstdata_obj,MSLDEMdata_obj,varargin)
            global msl_env_vars
            obj.cachedirpath = msl_env_vars.dirpath_cache;
            vr = ''; 
            % version 0: corresponds to original
            %         1: fixed rover nav
            if (rem(length(varargin),2)==1)
                error('Optional parameters should always go by pairs');
            else
                for i=1:2:(length(varargin)-1)
                    switch upper(varargin{i})
                        case {'VERSION','VER'}
                            vr = lower(varargin{i+1});
                        otherwise
                            error('Unrecognized option: %s',varargin{i});
                    end
                end
            end
            
            if isempty(vr)
                error('Please enter version with the option "VER" or "CACHE_VERSION".'); 
            else
                obj.version = vr;
            end
            
            obj.MASTCAMdata = mstdata_obj;
            obj.MSLDEMdata  = MSLDEMdata_obj;
            
            % setup cache file name.
            [basename_cache_com] = mastcam_create_basename_cache(mstdata_obj);
            obj.cachebasename_common = basename_cache_com;
            
            obj.rover_nav_version = mstdata_obj.ROVER_NAV.version;
            
        end
        
        % Master method performing all the projection
        function project(obj,varargin)
            % future implementation
        end
        
        %% dealing with cachefiles
        function loadCache(obj,varargin)
            obj.loadCache_projmastcam2MSLDEM(varargin{:});
            obj.loadCache_UFOVmask(varargin{:});
        end
        
        function loadCache_projmastcam2MSLDEM(obj,varargin)
            % this cache loader loads the file from UFOV cache data.
            fname_cache = sprintf('%s_projmastcam2MSLDEM_%s.mat',obj.cachebasename_common,obj.version);
            cachepath = joinPath(obj.cachedirpath,fname_cache);
            if exist(cachepath,'file')
                 [obj.mastcam_NEE,obj.mastcam_msldemc_ref,obj.mastcam_range,...
                     obj.mastcam_msldemc_nn,~,...
                     obj.mastcam_emi,obj.mastcam_surfplnc] = ...
                     mastcam_get_projmastcam2MSLDEM(obj.MASTCAMdata,...
                     obj.MSLDEMdata,obj,varargin{:},...
                     'Cache_Ver',obj.version);
            else
                error('No cache file exist');
            end
        end
        
        function loadCache_UFOVmask(obj,varargin)
            % this cache loader loads the file from UFOV cache data.
            basename_cache_ufovc = sprintf('%s_imUFOV_%s.mat',obj.cachebasename_common,obj.version);
            fpath_ufovc = joinPath(obj.cachedirpath,basename_cache_ufovc);
            if exist(fpath_ufovc,'file')
                crop_and_save_msldemc_imUFOVmask(obj,varargin{:},...
                    'Cache_Ver',obj.version);
            else
                error('No cache file exist');
            end
        end
        
        %% Projection component functions        
        function proj_MSLDEM2mastcam(obj,varargin)
            [obj.msldemc_imFOVmask,obj.msldemc_imFOVxy,obj.msldemc_imFOVhdr]...
                = mastcam_get_projMSLDEM2mastcam(obj.MSLDEMdata,obj.MASTCAMdata,varargin{:});
        end
        
        function proj_mastcam2MSLDEM(obj,varargin)
            [obj.mastcam_NEE,obj.mastcam_msldemc_ref,obj.mastcam_range,...
                obj.mastcam_msldemc_nn,obj.msldemc_imFOVmask_ctrnn,...
                obj.mastcam_emi,obj.mastcam_surfplnc] = ...
                mastcam_get_projmastcam2MSLDEM(obj.MASTCAMdata,...
                obj.MSLDEMdata,obj,varargin{:});
        end
        
        function getUFOVwMSLDEM(obj,varargin)
            % Optional input parameters can be 'BORDER_ASSESS_OPT',
            % 'CACHE_DIRPATH_CACHE','SAVE_FILE','Force'
            % 'LOAD_CACHE_IFEXIST', 'CACHE_VERSION'
            obj.msldemc_imUFOVmask = mastcam_getUFOVwMSLDEM_mexw(obj,varargin{:});
        end
        
        function crop_and_save_msldemc_imUFOVmask(obj,varargin)
             [obj.msldemc_imUFOVmask,obj.msldemc_imUFOVhdr,...
                 obj.msldemc_imUFOVxynn,obj.mastcam_msldemc_nn_UFOVmask]...
                 = mastcam_crop_msldemc_imUFOVmask(obj,varargin{:});
        end
        
        function get_FOVcone(obj,varargin)
            cmmdl_geo = transform_CAHVOR_MODEL_wROVER_NAV(obj.MASTCAMdata.CAM_MDL,obj.MASTCAMdata.ROVER_NAV);
            cmmdl_geo.get_image_plane_unit_vectors();
            L_im = obj.MASTCAMdata.L_im; S_im = obj.MASTCAMdata.S_im;
            [pmc_ul,pmc_ll,pmc_lr,pmc_ur] = cmmdl_geo.get_pmc_FOVvertex([0 S_im-1],[0 L_im-1]);
            
            Ane = cmmdl_geo.A(1:2) ./ sqrt(sum(cmmdl_geo.A(1:2).^2));
            pmc_ul_ne = pmc_ul(1:2) ./ sqrt(sum(pmc_ul(1:2).^2));
            pmc_ll_ne = pmc_ll(1:2) ./ sqrt(sum(pmc_ll(1:2).^2));
            pmc_lr_ne = pmc_lr(1:2) ./ sqrt(sum(pmc_lr(1:2).^2));
            pmc_ur_ne = pmc_ur(1:2) ./ sqrt(sum(pmc_ur(1:2).^2));
            if abs(sum(Ane*pmc_ul_ne)) > abs(sum(Ane*pmc_ll_ne))
                lvec = pmc_ll';
            else
                lvec = pmc_ul';
            end
            
            if abs(sum(Ane*pmc_lr_ne)) > abs(sum(Ane*pmc_ur_ne))
                rvec = pmc_ur';
            else
                rvec = pmc_lr';
            end
            
            mrange = obj.mastcam_range;
            mrange(isinf(mrange)) = nan;
            max_dst = max(mrange(:));
            
            cone_l = [cmmdl_geo.C;cmmdl_geo.C+lvec*max_dst];
            cone_r = [cmmdl_geo.C;cmmdl_geo.C+rvec*max_dst];
            obj.FOVcone.center = cmmdl_geo.C;
            obj.FOVcone.left   = cone_l;
            obj.FOVcone.right  = cone_r;
            
        end
        
        function plot_FOVcone(obj,ax,varargin)
            current_nextplot = ax.NextPlot;
            ax.NextPlot = 'add';
            plot(ax,obj.FOVcone.left(:,2),obj.FOVcone.left(:,1),'-');
            plot(ax,obj.FOVcone.right(:,2),obj.FOVcone.right(:,1),'-');
            plot(ax,obj.FOVcone.center(2),obj.FOVcone.center(1),'X');
            ax.NextPlot = current_nextplot;
            
        end
        
        function get_UFOVpolygons_quick(obj,varargin)
            % First surface plane and emission angnles are modified so that
            % emission anlges gets smaller than threshold values. Without
            % this, estimated resolution at some pixels may be so large.
            [mastcam_surfplnc_th,mastcam_emi_th] = mastcam_safeguard_planenormal(...
                obj.MASTCAMdata,obj.mastcam_emi,obj.mastcam_surfplnc,...
                obj.MASTCAMdata.ROVER_NAV);
            
            % Next, resolution in the H and V directions are calculated for
            % each pixels. The resolution is measured on the surface.
            [ifovH,ifovV] = mastcam_get_pixel_resolution(obj.MASTCAMdata,...
                obj.mastcam_range,mastcam_emi_th,mastcam_surfplnc_th,...
                obj.MASTCAMdata.ROVER_NAV);
            
            % Mastcam image is segmented 
            [impolid] = mastcam_get_FOVpolygons(obj.mastcam_NEE,2.5*ifovV,2.5*ifovH);
            
            % obtained segmentation is oraganized.
            [impolid_cleaned] = mastcam_clean_FOVpolygons(impolid,obj.mastcam_range);
            obj.qmastcam_imxypolygons.raster = impolid_cleaned;
            % Polygons are obtained. First polygons in the mastcam image
            [impolygoncell] = mastcam_get_FOVpolygons_Vectors(impolid_cleaned);
            obj.qmastcam_imxypolygons.vector = impolid_cleaned;
            % Then polygons in north-east-elevation domains.
            [obj.qNEEpolygons] = mastcam_get_NEEpolygons_from_impolygons(obj.mastcam_NEE,impolygoncell);
            
        end
        
        function plot_NEEpolygons(obj,ax,varargin)
            Np = length(obj.qNEEpolygons);
            current_nextplot = ax.NextPlot;
            ax.NextPlot = 'add';
            for i = 1:Np
                plot(ax,obj.qNEEpolygons{i}(:,2),obj.qNEEpolygons{i}(:,1),varargin{:});
            end
            ax.NextPlot = current_nextplot;
            
        end
        
        
        
        % function proj_MSLDEM2mastcam_old(obj)
        %     [MSLDEMprj] = proj_MSLDEM2mastcam_v2(obj.MSLDEMdata,obj.MASTCAMdata);
        %     l1 = MSLDEMprj.hdr_imxy.line_offset+1;
        %     lend = MSLDEMprj.hdr_imxy.line_offset+MSLDEMprj.hdr_imxy.lines;
        %     s1 = MSLDEMprj.hdr_imxy.sample_offset+1;
        %     send = MSLDEMprj.hdr_imxy.sample_offset+MSLDEMprj.hdr_imxy.samples;
        %     dem_imFOV_mask_crop = MSLDEMprj.imFOV_mask(l1:lend,s1:send);
        %     obj.msldemc_imFOVxy = MSLDEMprj.imxy;
        %     obj.msldemc_imFOVmask = dem_imFOV_mask_crop;
        %     obj.msldemc_imFOVhdr = MSLDEMprj.hdr_imxy;
        % end
        %         function crop_msldemc_imUFOVmask(obj)
%             valid_lines   = find(any(obj.msldemc_imUFOVmask',1));
%             lrnge         = [valid_lines(1), valid_lines(end)];
%             len_vl        = lrnge(2)-lrnge(1)+1;
%             valid_samples = find(any(obj.msldemc_imUFOVmask,1));
%             srnge         = [valid_samples(1), valid_samples(end)];
%             len_vs        = srnge(2)-srnge(1)+1;
% 
%             l1   = obj.msldemc_imFOVhdr.line_offset+lrnge(1);
%             lend = obj.msldemc_imFOVhdr.line_offset+lrnge(2);
%             s1   = obj.msldemc_imFOVhdr.sample_offset+srnge(1);
%             send = obj.msldemc_imFOVhdr.sample_offset+srnge(2);
%             msldemc_northing = obj.MSLDEMdata.hdr.y(l1:lend);
%             msldemc_easting  = obj.MSLDEMdata.hdr.x(s1:send);
% 
%             obj.msldemc_imUFOVmask = obj.msldemc_imUFOVmask(lrnge(1):lrnge(2),srnge(1):srnge(2));
% 
%             obj.msldemc_imUFOVhdr = [];
%             obj.msldemc_imUFOVhdr.lines   = len_vl;
%             obj.msldemc_imUFOVhdr.samples = len_vs;
%             obj.msldemc_imUFOVhdr.line_offset   = l1-1;
%             obj.msldemc_imUFOVhdr.sample_offset = s1-1;
%             obj.msldemc_imUFOVhdr.y = msldemc_northing;
%             obj.msldemc_imUFOVhdr.x = msldemc_easting;
% 
%             obj.msldemc_imUFOVxy = obj.msldemc_imFOVxy(lrnge(1):lrnge(2),srnge(1):srnge(2),:);
%             mm = (obj.msldemc_imUFOVmask==0);
%             for i=1:2
%                 msldemc_imUFOVtmp = obj.msldemc_imUFOVxy(:,:,i);
%                 msldemc_imUFOVtmp(mm) = nan;
%                 obj.msldemc_imUFOVxy(:,:,i) = msldemc_imUFOVtmp;
%             end
%             
%             obj.msldemc_imUFOVxynn = round(obj.msldemc_imUFOVxy+1);
%             
%             L_im = obj.MASTCAMdata.L_im; S_im = obj.MASTCAMdata.S_im;
%             
%             msldemc_imUFOVtmp = obj.msldemc_imUFOVxynn(:,:,1);
%             msldemc_imUFOVtmp(msldemc_imUFOVtmp<1) = 1;
%             msldemc_imUFOVtmp(isnan(msldemc_imUFOVtmp)) = -1;
%             msldemc_imUFOVtmp(msldemc_imUFOVtmp>obj.MASTCAMdata.S_im) = S_im;
%             obj.msldemc_imUFOVxynn(:,:,1) = msldemc_imUFOVtmp;
%             msldemc_imUFOVtmp = obj.msldemc_imUFOVxynn(:,:,2);
%             msldemc_imUFOVtmp(msldemc_imUFOVtmp<1) = 1;
%             msldemc_imUFOVtmp(isnan(msldemc_imUFOVtmp)) = -1;
%             msldemc_imUFOVtmp(msldemc_imUFOVtmp>obj.MASTCAMdata.L_im) = L_im;
%             obj.msldemc_imUFOVxynn(:,:,2) = msldemc_imUFOVtmp;
%             obj.msldemc_imUFOVxynn = int16(obj.msldemc_imUFOVxynn);
%             
%             obj.mastcam_msldemc_nn_UFOVmask = obj.mastcam_msldemc_nn;
%             obj.mastcam_msldemc_nn_UFOVmask(:,:,1) = obj.mastcam_msldemc_nn_UFOVmask(:,:,1)-int32(srnge(1)-1);
%             obj.mastcam_msldemc_nn_UFOVmask(:,:,2) = obj.mastcam_msldemc_nn_UFOVmask(:,:,2)-int32(lrnge(1)-1);
%             
%         end
        %% Interfaces for resolving pixel matching
        function [x_east,y_north] = get_NE_from_msldemc_imUFOVxy(obj,x_dem,y_dem)
            x_east = obj.msldemc_imUFOVhdr.x(x_dem);
            y_north = obj.msldemc_imUFOVhdr.y(y_dem);
        end
        
        function [x_dem,y_dem] = get_imUFOVc_xy_from_NE(obj,x_east,y_north)
            [x_dem] = round((x_east-obj.msldemc_imUFOVhdr.x(1))/obj.MSLDEMdata.hdr.map_info.dx + 1);
            [y_dem] = round((obj.msldemc_imUFOVhdr.y(1)-y_north)/obj.MSLDEMdata.hdr.map_info.dy + 1);
            if x_dem < 1 || x_dem > obj.msldemc_imUFOVhdr.samples || y_dem < 1 || y_dem > obj.msldemc_imUFOVhdr.lines
                x_dem = nan; y_dem = nan;
            end
        end
        
        function [x_mst,y_mst] = get_mastcam_imxy_from_imUFOVcxy(obj,x_dem,y_dem)
            x_mst = obj.msldemc_imUFOVxynn(y_dem,x_dem,1);
            y_mst = obj.msldemc_imUFOVxynn(y_dem,x_dem,2);
        end
        
        function [x_dem,y_dem] = get_msldemc_imUFOVxy_from_mastcam_imxy(obj,x_mst,y_mst)
            x_dem = obj.mastcam_msldemc_nn_UFOVmask(y_mst,x_mst,1);
            y_dem = obj.mastcam_msldemc_nn_UFOVmask(y_mst,x_mst,2);
        end
        
        function [x_east,y_north] = get_NorthEast_from_mastcam_imxy(obj,x_mst,y_mst)
            x_east = obj.mastcam_NEE(y_mst,x_mst,2);
            y_north = obj.mastcam_NEE(y_mst,x_mst,1);
        end
        
        function [x_mst,y_mst,x_dem,y_dem,mask_mastcam,mask_msldemc_UFOVmask]...
                = get_SelectMask_fromNorthEast(obj,x_east,y_north)
            % [x_mst,y_mst,x_dem,y_dem,mask_mastcam,mask_msldemc_UFOVmask]...
            %     = obj.get_SelectMask_fromNorthEast(x_east,y_north)
            % get select mask give (east, north) input
            %  INPUTS
            %    x_east, y_north: meters, easting and northing. 
            %  OUTPUTS
            %    x_mst, y_mst : corresponding central mastcam image coordinates
            %    x_dem, y_dem : nearest neighbor of the image pixel 
            %       coordinate in msldemc_imUFOVmask
            %    mask_mastcam : L_im x S_im boolean image, true if the
            %    mastcam image pixel corresponds to (x_dem,y_dem) in
            %    msldemc_imUFOVmask.
            %    mask_msldemc_UFOVmask : boolean image, corresponding
            %    pixels (x_mst,y_mst) in the mastcam image.
            [x_dem,y_dem] = obj.get_imUFOVc_xy_from_NE(x_east,y_north);
            if ~isnan(x_dem) && ~isnan(y_dem) && (obj.msldemc_imUFOVmask(y_dem,x_dem)>0)
                [x_mst,y_mst] = obj.get_mastcam_imxy_from_imUFOVcxy(x_dem,y_dem);
                
                % MASKs (forward and backward)
                mask_msldemc_UFOVmask = false(obj.msldemc_imUFOVhdr.lines,obj.msldemc_imUFOVhdr.samples);
                mask_msldemc_UFOVmask(y_dem,x_dem) = true;
                idxes = obj.mapper_msldemc2mastcam_cell{obj.mapper_msldemc2mastcam_mat(y_dem,x_dem)};
                mask_mastcam = false(obj.MASTCAMdata.L_im,obj.MASTCAMdata.S_im);
                for i=1:size(idxes,2)
                    mask_mastcam(idxes(2,i),idxes(1,i)) = true;
                    idxes_j = obj.mapper_mastcam2msldemc{idxes(2,i),idxes(1,i)};
                    for j=1:size(idxes_j,2)
                        mask_msldemc_UFOVmask(idxes_j(2,j),idxes_j(1,j)) = true;
                    end
                end
                % MASK for MSLDEM image
                
                
                % [mask_mastcam] = obj.create_mask_mastcam_wMSLDEMcoord(...
                %     x_dem,y_dem,x_mst,y_mst);
                % [mask_msldemc_UFOVmask] = obj.create_mask_msldem_wMSLDEMcoord(...
                %     x_dem,y_dem,x_mst,y_mst);
            else
                x_mst = nan; y_mst = nan;
                mask_mastcam = []; mask_msldemc_UFOVmask = [];
            end
        end
        
        function [x_east,y_north,x_dem,y_dem,mask_mastcam,mask_msldemc_UFOVmask]...
                = create_SelectMask_fromMASTCAM(obj,x_mst,y_mst)
            
            % MASKs (forward and backward)
            mask_mastcam = false(obj.MASTCAMdata.L_im,obj.MASTCAMdata.S_im);
            mask_mastcam(y_mst,x_mst) = true;
            idxes = obj.mapper_mastcam2msldemc{y_mst,x_mst};
            mask_msldemc_UFOVmask = false(obj.msldemc_imUFOVhdr.lines,obj.msldemc_imUFOVhdr.samples);

            for i=1:size(idxes,2)
                mask_msldemc_UFOVmask(idxes(2,i),idxes(1,i)) = true;
                idxes_j = obj.mapper_msldemc2mastcam_cell{obj.mapper_msldemc2mastcam_mat(idxes(2,i),idxes(1,i))};
                for j=1:size(idxes_j,2)
                    mask_mastcam(idxes_j(2,j),idxes_j(1,j)) = true;
                end
            end
            
            % 
            [x_east,y_north] = obj.get_NorthEast_from_mastcam_imxy(x_mst,y_mst);
            [x_dem,y_dem] = obj.get_msldemc_imUFOVxy_from_mastcam_imxy(x_mst,y_mst);
            % if ~isnan(x_dem) && ~isnan(y_dem) && (obj.msldemc_imUFOVmask(y_dem,x_dem)>0)
            %     [x_mst,y_mst] = obj.get_mastcam_imxy_from_imUFOVcxy(x_dem,y_dem);
            %     
            %     % MASK for MASTCAM image using MSLDEM
            %     [mask_mastcam] = obj.create_mask_mastcam_wMSLDEMcoord(...
            %         x_dem,y_dem,x_mst,y_mst);
            %     [mask_msldemc_UFOVmask] = obj.create_mask_msldem_wMSLDEMcoord(...
            %         x_dem,y_dem,x_mst,y_mst);
            % else
            %     mask_mastcam = []; mask_msldemc_UFOVmask = [];
            % end
        end
        
        function [mask_mastcam,mask_msldemc_UFOVmask]...
                = create_SelectMask_from_msldemc_mask(obj,msldemc_mask)
            
            % MASKs (forward and backward)
            mask_msldemc_UFOVmask = msldemc_mask;
            
            mask_mastcam = false(obj.MASTCAMdata.L_im,obj.MASTCAMdata.S_im);
            
            [ufov_row,ufov_col] = find(msldemc_mask);
            for i=1:length(ufov_row)
                ufov_rowi = ufov_row(i); ufov_coli = ufov_col(i);
                idxes = obj.mapper_msldemc2mastcam_cell{obj.mapper_msldemc2mastcam_mat(ufov_rowi,ufov_coli)};
                for ii=1:size(idxes,2)
                    mask_mastcam(idxes(2,ii),idxes(1,ii)) = true;
                    % next is back projection
                    idxes_j = obj.mapper_mastcam2msldemc{idxes(2,ii),idxes(1,ii)};
                    for j=1:size(idxes_j,2)
                        mask_msldemc_UFOVmask(idxes_j(2,j),idxes_j(1,j)) = true;
                    end
                end
            end
            
        end
        
%         function [mask_mastcam] = create_mask_mastcam_wMSLDEMcoord(obj,x_dem,y_dem,x_mst,y_mst)
%             % MASK for MASTCAM image using MSLDEM
%             % select all the pixels of MASTCAM image whose nearest neighbor 
%             % is (x_dem,y_dem) in the MSLDEM image space.
%             % You should get nozeros when the selected MSLDEM image pixel
%             % covers multiple pixels of MASTCAM image.
%             % --- SINGLE MSLDEM pixel -------------------------------------
%             mask_mastcam = and(...
%                 obj.mastcam_msldemc_nn_UFOVmask(:,:,1)==x_dem,...
%                 obj.mastcam_msldemc_nn_UFOVmask(:,:,2)==y_dem);
%             % The nearest MASTCAM pixel from the picked MSLDEM image pixel
%             % is also selected.
%             mask_mastcam(y_mst,x_mst) = true;
%         end
%         
%         function [mask_msldemc_UFOVmask] = create_mask_msldem_wMSLDEMcoord(obj,x_dem,y_dem,x_mst,y_mst)
%             % MASK for MSLDEM image
%             % select all the pixels of MSLDEM image whose nearest neighbor 
%             % is (x_mst,y_mst) in the mastcam image space.
%             % You should get nozeros when the selected MASTCAM image pixel
%             % covers multiple pixels of MSLDEM image.
%             % This mask might cover multiple MSLDEM image pixels if the
%             % pointed MASTCAM pixel covers muttiple MSLDEM image pixels.
%             % mask_msldemc_UFOVmask = and(obj.msldemc_imUFOVxynn(:,:,1)==x_mst,...
%             %     obj.msldemc_imUFOVxynn(:,:,2)==y_mst);
%             mask_msldemc_UFOVmask = false(obj.msldemc_imUFOVhdr.lines,obj.msldemc_imUFOVhdr.samples);
%             for i=1:size(obj.mapper_mastcam2msldemc{y_mst,x_mst},2)
%                 tmp = obj.mapper_mastcam2msldemc{y_mst,x_mst};
%                 mask_msldemc_UFOVmask(tmp(2,i),tmp(1,i)) = true;
%             end
%             mask_msldemc_UFOVmask(y_dem,x_dem) = true;
%         end
        
        function [msldemcUFOV_mask,mask_mastcam] = get_rangeUFOVmask(obj,xrnge_east,yrnge_north)
            xUFOV_within =  and(obj.msldemc_imUFOVhdr.x > xrnge_east(1),obj.msldemc_imUFOVhdr.x < xrnge_east(2));
            yUFOV_within =  and(obj.msldemc_imUFOVhdr.y > yrnge_north(1),obj.msldemc_imUFOVhdr.y < yrnge_north(2));
            msldemcUFOV_mask = false(obj.msldemc_imUFOVhdr.lines,obj.msldemc_imUFOVhdr.samples);
            msldemcUFOV_mask(yUFOV_within,xUFOV_within) = true;
            
            msldemcUFOV_mask = and(msldemcUFOV_mask,obj.msldemc_imUFOVmask);
            
            [ufov_row,ufov_col] = find(msldemcUFOV_mask);
            
            mask_mastcam = false(obj.MASTCAMdata.L_im,obj.MASTCAMdata.S_im);
            tic;
            for i=1:length(ufov_row)
                ufov_rowi = ufov_row(i); ufov_coli = ufov_col(i);
                [x_mst,y_mst] = obj.get_mastcam_imxy_from_imUFOVcxy(ufov_coli,ufov_rowi);
                [mask_mastcam_i] = obj.create_mask_mastcam_wMSLDEMcoord(ufov_coli,ufov_rowi,x_mst,y_mst);
                [msldemcUFOV_mask_i] = obj.create_mask_msldem_wMSLDEMcoord(ufov_coli,ufov_rowi,x_mst,y_mst);
                mask_mastcam = or(mask_mastcam,mask_mastcam_i);
                tic; msldemcUFOV_mask = or(msldemcUFOV_mask,msldemcUFOV_mask_i); toc;
            end
            toc;
        end
        
        %%
        function [msldemc_img,msldemc_img_hdr] = MSLDEMread(obj)
            mrgn = 100;
            smpl_offset = max(1,obj.msldemc_imUFOVhdr.sample_offset-mrgn);
            smpls = min(obj.MSLDEMdata.hdr.samples-smpl_offset,...
                obj.msldemc_imUFOVhdr.samples+2*mrgn);
            ln_offset = max(1,obj.msldemc_imUFOVhdr.line_offset-mrgn);
            lns = min(obj.MSLDEMdata.hdr.lines-ln_offset,...
                obj.msldemc_imUFOVhdr.lines+2*mrgn);
            
            tic; msldemc_img = msldem_lazyenvireadRect(obj.MSLDEMdata,...
                smpl_offset,ln_offset,smpls,lns,'precision','single'); toc;
            
            msldemc_img_hdr = [];
            msldemc_img_hdr.sample_offset = smpl_offset;
            msldemc_img_hdr.line_offset   = ln_offset;
            msldemc_img_hdr.samples = smpls;
            msldemc_img_hdr.lines   = lns;
            
        end
    end
end

