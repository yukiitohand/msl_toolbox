classdef MASTCAMCameraProjectionMSLDEM_v2 < handle
    %MASTCAMCameraProjectionMSLDEM_v2
    %  class handling camera projection.
    %  Properties
    %    MASTCAMdata: MASTCAMdata or MASTCAMgroup_eye class object
    %    MSLDEMdata: MSLGaleDEMMosaic_v3 class obj, 
    %      MSLDEMdata at
    %      https://astrogeology.usgs.gov/search/map/Mars/MarsScienceLaboratory/
    %      Mosaics/MSL_Gale_DEM_Mosaic_10m
    %      or its variants
    %
    %    
    %     
    
    properties
        MASTCAMdata
        % MASTCAMdata or MASTCAMgroup_wProcCode class object
        
        MSLDEMdata
        % MSLGaleDEMMosaic_v3 class or its subclass object
        
        msldemc_imFOV
        % MSLDEMCimFOV class object
        %  Properties
        %    mask
        %    maskr13
        %    pxlctrnn
        
        msldemc_imUFOV
        % MSLDEMCimUFOV class object
        %  Properties
        %    mask
        %    xynn
        
        mastcam_proj
        % MASTCAMgroup_projection class object
        %  Properties
        %   XYZ
        %   NE
        %   latlon
        %   zenith
        %   range
        %   surfplnc
        %   emiang
        %   nn_msldem
        %   ref_msldem
        
        mapper
        % MASTCAM_MSLDEMC_Mapper class object
        
        cachedirpath
        cachebasename_common
        version
        rover_nav_version
        FOVcone
        qNEEpolygons          % quick NEE polygons
        qmastcam_imxypolygons % quick imxy polygons
        
        identifier;
        
        
        
    end
    
    methods
        function obj = MASTCAMCameraProjectionMSLDEM_v2(mstdata_obj,MSLDEMdata_obj,varargin)
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
                            validateattributes(vr,{'char'},{}, ...
                                mfilename,'CACHE_VERSION');
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
            
            %
            % Input validation
            validateattributes(mstdata_obj, ...
                {'MASTCAMdata','MASTCAMgroup_wProcCode'},{}, ...
                mfilename,'mstdata_obj');
            validateattributes(MSLDEMdata_obj,{'MSLGaleDEMMosaic_v3'}, ...
                {},mfilename,'MSLDEMdata_obj');
            
            if isempty(mstdata_obj.PRODUCT_ID)
                error('No data is found');
            end
            
            obj.MASTCAMdata = mstdata_obj;
            obj.MSLDEMdata  = MSLDEMdata_obj;
            
            % setup cache file name.
            [basename_cache_com] = mastcam_create_basename_cache(mstdata_obj);
            obj.cachebasename_common = basename_cache_com;
            
            obj.rover_nav_version = mstdata_obj.ROVER_NAV.version;
            
            obj.identifier = basename_cache_com;
            
        end
        
        % Master method performing all the projection
        function project(obj,varargin)
            % future implementation
        end
        
        
        %% Projection component functions        
        function proj_MSLDEM2mastcam(obj,varargin)
            if isempty(obj.msldemc_imFOV)
                obj.msldemc_imFOV = MSLDEMCimFOV();
            end
            obj.msldemc_imFOV.mask = mastcam_get_msldemc_imFOVmask( ...
                obj.MSLDEMdata,obj.MASTCAMdata,varargin{:});
        end
        
        function proj_mastcam2MSLDEM(obj,varargin)
            if isempty(obj.msldemc_imFOV)
                obj.msldemc_imFOV = MSLDEMCimFOV();
            end
            [obj.mastcam_proj,obj.msldemc_imFOV.pxlctrnn] =            ...
                mastcam_get_projmastcam2MSLDEM_v2(obj.MASTCAMdata,     ...
                obj.MSLDEMdata,obj.msldemc_imFOV.mask,varargin{:});
        end
        
        function msldemc_imFOVmask_eval_13(obj,varargin)
            if isempty(obj.msldemc_imFOV)
                obj.msldemc_imFOV = MSLDEMCimFOV();
            end
            [obj.msldemc_imFOV.maskr13] = mastcam_get_msldemc_imFOVmask_13refined( ...
                obj.MASTCAMdata,obj.msldemc_imFOV.mask,obj.mastcam_proj.ref_msldem,varargin{:});
            % [obj.msldemc_imFOVmask] = mastcam_msldemc_imFOVmask_eval_13(obj.msldemc_imFOVmask,obj.mastcam_msldemc_ref);
        end
        
        function getUFOVwMSLDEM(obj,varargin)
            % Optional input parameters can be 'BORDER_ASSESS_OPT',
            % 'CACHE_DIRPATH_CACHE','SAVE_FILE','Force'
            % 'LOAD_CACHE_IFEXIST', 'CACHE_VERSION'
            if isempty(obj.msldemc_imUFOV)
                obj.msldemc_imUFOV = MSLDEMCimUFOV();
            end
            obj.msldemc_imUFOV.mask = mastcam_get_imUFOVmask(...
                obj.MASTCAMdata,obj.MSLDEMdata, ...
                obj.msldemc_imFOV.maskr13,obj.msldemc_imFOV.pxlctrnn,varargin{:});
        end
        
        function get_mapper(obj,varargin)
            if isempty(obj.msldemc_imUFOV)
                obj.msldemc_imUFOV = MSLDEMCimUFOV();
            end
            if isempty(obj.mastcam_proj)
                obj.mastcam_proj = MASTCAMgroup_projection();
            end
            obj.msldemc_imUFOV.xynn = mastcam_get_imUFOVxynn(obj.MASTCAMdata,obj.MSLDEMdata,obj.msldemc_imUFOV.mask,varargin{:});
            obj.mapper = mastcam_get_msldemc_mapper(obj.MASTCAMdata,obj.MSLDEMdata,obj.msldemc_imUFOV.xynn,obj.mastcam_proj.nn_msldem,varargin{:});
        end
        
        %%
        
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
        
        %% Interfaces for resolving pixel matching
        function [x_east,y_north] = get_NE_from_msldemc_imUFOVxy(obj,x_dem,y_dem)
            x_east = obj.msldemc_imUFOVhdr.x(x_dem);
            y_north = obj.msldemc_imUFOVhdr.y(y_dem);
        end
        
        function [x_dem,y_dem] = get_imUFOVc_xy_from_NE(obj,x_east,y_north)
            [x_dem] = round((x_east-obj.msldemc_imUFOVhdr.x(1))/obj.MSLDEMdata.proj_info.map_scale_x + 1);
            [y_dem] = round((obj.msldemc_imUFOVhdr.y(1)-y_north)/obj.MSLDEMdata.proj_info.map_scale_y + 1);
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
            
            xrange = [smpl_offset+1 smpl_offset+smpls];
            yrange = [ln_offset+1 ln_offset+lns];
            
            tic; msldemc_img = obj.MSLDEMdata.get_subimage_wPixelRange(...
                xrange,yrange,'precision','single'); toc;
            
            msldemc_img_hdr = [];
            msldemc_img_hdr.sample_offset = smpl_offset;
            msldemc_img_hdr.line_offset   = ln_offset;
            msldemc_img_hdr.samples = smpls;
            msldemc_img_hdr.lines   = lns;
            
        end
    end
end

