classdef MASTCAMprojectionView_wMSLDEM < handle
    %MASTCAMprojectionView_wMSLDEM
    %   Viewer for MASTCAM images with projection onto MSLDEM
    
    properties
        MSTview
        ISV_MSLDEM
        MSTproj
        mstdata
        objHSIview
        MSTview_ISVImage_range
        isvimage_NEE_SelectMask
        isvimage_MASTCAM_SelectMask
    end
    
    methods
        function obj = MASTCAMprojectionView_wMSLDEM(MSTproj,varargin)
            %Constructor
            if ~isempty(varargin)
                for i=1:2:length(varargin)
                    switch upper(varargin{i})
                        % case {'IMAGE_NAMES','IMAGE_TITLES'}
                        %     obj.image_names = varargin{i+1};
                        % case 'CLIM'
                        %     image_clims = varargin{i+1};
                        case 'HSIVIEW'
                            obj.objHSIview = varargin{i+1};
                        otherwise
                            error('Unrecognized option: %s', varargin{i});
                    end
                end
            end
            %   
            obj.MSTproj = MSTproj;
            obj.MSTview = MASTCAMview(MSTproj.MASTCAMdata);
            obj.MSTview_ISVImage_range = ...
                obj.MSTview.obj_HSIview.obj_ISV.add_layer(...
                MSTproj.mastcam_range,'name','Range');
            obj.MSTview.obj_HSIview.obj_ISV.custom_image_cursor_fcn = @obj.ISV_MASTCAM_BtnDwnFcn;
            obj.MSTview.obj_HSIview.obj_ISV.custom_windowkeypress_fcn = @obj.ISV_MASTCAM_WindowKeyPressFcn;
            
            if isempty(obj.objHSIview)
                obj.ISV_MSLDEM = ImageStackView({},'Ydir','normal',...
                    'XY_COORDINATE_SYSTEM','NorthEast',...
                    'IMAGE_CURSOR_FCN',@obj.ISV_MSLDEM_BtnDwnFcn,...
                    'IMAGE_WINDOWKEYPRESS_FCN',@obj.ISV_MSLDEM_WindowKeyPressFcn);
                obj.init_ISV_MSLDEM();
            else
                obj.ISV_MSLDEM = obj.objHSIview.obj_ISV;
                obj.ISV_MSLDEM.custom_image_cursor_fcn = @obj.ISV_MSLDEM_BtnDwnFcn;
                obj.ISV_MSLDEM.custom_windowkeypress_fcn = @obj.ISV_MSLDEM_WindowKeyPressFcn;
                obj.init_ISV_MSLDEM();
            end
        end
        
        function init_ISV_MSLDEM(obj)
            % get MSLDEM image
            % set margin pixels for the UFOVmask
            mrgn = 1000;
            smpl_offset = max(1,obj.MSTproj.msldemc_imUFOVhdr.sample_offset-mrgn);
            smpls = min(obj.MSTproj.MSLDEMdata.hdr.samples-smpl_offset,...
                obj.MSTproj.msldemc_imUFOVhdr.samples+2*mrgn);
            ln_offset = max(1,obj.MSTproj.msldemc_imUFOVhdr.line_offset-mrgn);
            lns = min(obj.MSTproj.MSLDEMdata.hdr.lines-ln_offset,...
                obj.MSTproj.msldemc_imUFOVhdr.lines+2*mrgn);
            
            tic; msldemc_img = msldem_lazyenvireadRect(obj.MSTproj.MSLDEMdata,...
                smpl_offset,ln_offset,smpls,lns,'precision','single'); toc;           
            
            obj.ISV_MSLDEM.add_layer(...
                obj.MSTproj.MSLDEMdata.hdr.x([smpl_offset+1,smpl_offset+smpls]),...
                obj.MSTproj.MSLDEMdata.hdr.y([ln_offset+1,ln_offset+lns]),...
                msldemc_img,'name','MSLDEM');
            obj.ISV_MSLDEM.add_layer(...
                obj.MSTproj.msldemc_imUFOVhdr.x([1,end]),...
                obj.MSTproj.msldemc_imUFOVhdr.y([1,end]),...
                obj.MSTproj.msldemc_imUFOVmask,'name','UFOV mask',...
                'alphadata',obj.MSTproj.msldemc_imUFOVmask>0);
            
            obj.ISV_MSLDEM.Update_axim_aspectR();
            obj.ISV_MSLDEM.Restore_ImageAxes2LimHome();

        end
        
        function init_ISV_MSLDEMquick(obj)
            % get MSLDEM image
            % set margin pixels for the UFOVmask
            mrgn = 1000;
            north_min = min(min(obj.MSTproj.mastcam_NEE(:,:,1),[],'all'),obj.MSTproj.FOVcone.center(1));
            north_max = max(max(obj.MSTproj.mastcam_NEE(:,:,1),[],'all'),obj.MSTproj.FOVcone.center(1));
            
            y_strt = round(abs(obj.MSTproj.MSLDEMdata.hdr.map_info.mapy - north_max) ./ obj.MSTproj.MSLDEMdata.hdr.map_info.dy);
            y_end =  round(abs(obj.MSTproj.MSLDEMdata.hdr.map_info.mapy - north_min) ./ obj.MSTproj.MSLDEMdata.hdr.map_info.dy);
            
            east_min = min(min(obj.MSTproj.mastcam_NEE(:,:,2),[],'all'),obj.MSTproj.FOVcone.center(2));
            east_max = max(max(obj.MSTproj.mastcam_NEE(:,:,2),[],'all'),obj.MSTproj.FOVcone.center(2));
            
            x_strt = round(abs(east_min - obj.MSTproj.MSLDEMdata.hdr.map_info.mapx) ./ obj.MSTproj.MSLDEMdata.hdr.map_info.dx);
            x_end =  round(abs(east_max - obj.MSTproj.MSLDEMdata.hdr.map_info.mapx) ./ obj.MSTproj.MSLDEMdata.hdr.map_info.dx);
            
            
            
            smpl_offset = max(1,x_strt-mrgn);
            smpls = min(x_end-x_strt+1+mrgn*2, obj.MSTproj.MSLDEMdata.hdr.samples-smpl_offset);
            ln_offset = max(1,y_strt-mrgn);
            lns = min(y_end-y_strt+1+2*mrgn, obj.MSTproj.MSLDEMdata.hdr.lines-ln_offset);
            
            tic; msldemc_img = msldem_lazyenvireadRect(obj.MSTproj.MSLDEMdata,...
                smpl_offset,ln_offset,smpls,lns,'precision','single'); toc;           
            
            obj.ISV_MSLDEM.add_layer(...
                obj.MSTproj.MSLDEMdata.hdr.x([smpl_offset+1,smpl_offset+smpls]),...
                obj.MSTproj.MSLDEMdata.hdr.y([ln_offset+1,ln_offset+lns]),...
                msldemc_img,'name','MSLDEM');
            
            ax = obj.ISV_MSLDEM.ax_plot;
            % ax.NextPlot = 'add';
            obj.MSTproj.plot_FOVcone(ax);
            obj.MSTproj.plot_NEEpolygons(ax);
            
            
            obj.ISV_MSLDEM.Update_axim_aspectR();
            obj.ISV_MSLDEM.Restore_ImageAxes2LimHome();

        end
        
        function [im_mastcam] = get_rgb_MASTCAM(obj)
            if isprop(obj.mstdata,'E')
                obj.mstdata.E.DRLX.readimg();
                im_mastcam = obj.mstdata.E.DRLX.img;
            elseif isprop(obj.mstdata,'D')
                obj.mstdata.D.DRLX.readimg();
                im_mastcam = obj.mstdata.D(1).DRLX.img;
            end
        end
        
        %%
        % =================================================================
        % MSLDEM button down callback function
        % =================================================================
        function [out] = ISV_MSLDEM_BtnDwnFcn(obj,hObject,eventData)
            if ~isempty(obj.objHSIview)
                [out] = obj.objHSIview.image_BtnDwnFcn_HSIview(hObject,eventData);
            else
                [out] = obj.ISV_MSLDEM.image_BtnDwnFcn(hObject,eventData);
            end
            cursor_obj = out.cursor_obj;
            if isfield(cursor_obj.UserData,'MSTprojViewPlotObj')
            else
                cursor_obj.UserData.MSTprojViewPlotObj = MSTprojViewPlot();
            end
            
            % Obtain coordinate in the MSLDEM image pixel domain.
            
            x_east = out.cursor_obj.X; y_north = out.cursor_obj.Y;
            % =============================================================
            % First get the selected point/area in the MASTCAM image at 
            % MSLDEM resolution.
            [x_mst,y_mst,x_dem,y_dem,mask_mastcam,mask_msldemc_UFOVmask]...
                = obj.MSTproj.get_SelectMask_fromNorthEast(x_east,y_north);

            % --- SINGLE HSI pixel ----------------------------------------
            % Evaluate HSI
            obj.ISV_MSLDEM_BtnDwnFcn_HSI(cursor_obj,x_east,y_north,...
                mask_msldemc_UFOVmask);
            
            % Show mask at MSLDEM resolution.
            % This is added after HSI image mask as we expect this would be
            % smaller than HSI footprints.
            if ~isempty(mask_mastcam)
                % nextplot_cur = obj.isvimage_MASTCAM_SelectMask.ax.NextPlot;
                % obj.isvimage_MASTCAM_SelectMask.ax.NextPlot = 'add';
                imobj = obj.update_isvimage_MASTCAM_SelectMask(mask_mastcam,[0.8,0.4,0]);
                cursor_obj.UserData.MSTprojViewPlotObj.imObj_msldemc_UFOVmask_msldemresol = imobj;
                % obj.isvimage_MASTCAM_SelectMask.ax.NextPlot = nextplot_cur;
                
                lnObj_spc_mastcam = obj.plot_mastcam(obj.MSTview,mask_mastcam);
                cursor_obj.UserData.MSTprojViewPlotObj.lnObj_spc_mastcam_msldemresol = lnObj_spc_mastcam;
            end
            if ~isempty(mask_msldemc_UFOVmask)
                % nextplot_cur = obj.isvimage_NEE_SelectMask.ax.NextPlot;
                % obj.isvimage_NEE_SelectMask.ax.NextPlot = 'add';
                imobj = obj.update_isvimage_NEE_SelectMask(mask_msldemc_UFOVmask,[0.8,0.4,0]);
                cursor_obj.UserData.MSTprojViewPlotObj.imObj_mastcam_msldemresol = imobj;
                % obj.isvimage_NEE_SelectMask.ax.NextPlot = nextplot_cur;
            end
            
        end
        
        function [out] = ISV_MSLDEM_WindowKeyPressFcn(obj,figobj,eventData)
            [out] = obj.ISV_MSLDEM.ISVWindowKeyPressFcn(figobj,eventData);
            % future implementation
        end

        function [] = ISV_MSLDEM_BtnDwnFcn_HSI(obj,cursor_obj,x_east,y_north,mask_msldemc_UFOVmask)
            for i=1:obj.objHSIview.nhsi
                hsielem = obj.objHSIview.hsiar(i);
                % get range of the selected hsi pixel.
                [mask_hsi,mask_msldemc_UFOVmask_hsiresol,...
                    mask_mastcam_hsiresol] = obj.get_hsi_mask(...
                    hsielem,x_east,y_north,mask_msldemc_UFOVmask);

                [imobj] = obj.update_isvimage_NEE_SelectMask(mask_msldemc_UFOVmask_hsiresol,[0,0,0.8]);
                cursor_obj.UserData.MSTprojViewPlotObj.imObj_msldemc_UFOVmask_hsiresol(i) = imobj;
                lnObj_spc_mastcam = obj.plot_mastcam(obj.MSTview,mask_mastcam_hsiresol);
                cursor_obj.UserData.MSTprojViewPlotObj.lnObj_spc_mastcam_hsiresol(i) = lnObj_spc_mastcam;
                
                [imobj] = obj.update_isvimage_MASTCAM_SelectMask(mask_mastcam_hsiresol,[0,0,0.8]);
                cursor_obj.UserData.MSTprojViewPlotObj.imObj_mastcam_hsiresol(i) = imobj;
                
            end
        end
        
        %%
        % =================================================================
        % MASTCAM window button down callback function
        % =================================================================
        function [out] = ISV_MASTCAM_BtnDwnFcn(obj,hObject,eventData)
            [out] = obj.MSTview.obj_HSIview.image_BtnDwnFcn_HSIview(hObject,eventData);
            cursor_obj = out.cursor_obj;
            if isfield(cursor_obj.UserData,'MSTprojViewPlotObj')
            else
                cursor_obj.UserData.MSTprojViewPlotObj = MSTprojViewPlot();
            end
            x_mst = out.cursor_obj.X; y_mst = out.cursor_obj.Y;
            
            % if cursor is pointing the sky, then do not any further
            % processing
            if isinf(obj.MSTproj.mastcam_range(y_mst,x_mst))
            else
                % First draw the selected point in the MASTCAM image at 
                % MSLDEM resolution.
                % First get the selected point/area in the MSLDEM image at 
                % MSLDEM resolution.
                tic; [x_east,y_north,x_dem,y_dem,mask_mastcam,mask_msldemc_UFOVmask]...
                    = obj.MSTproj.create_SelectMask_fromMASTCAM(x_mst,y_mst); toc;
                
                
                % --- SINGLE HSI pixel ----------------------------------------
                obj.ISV_MSLDEM_BtnDwnFcn_HSI(cursor_obj,x_east,y_north,mask_msldemc_UFOVmask);
                
                % Show mask at MSLDEM resolution.
                % This is added after HSI image mask as we expect this would be
                % smaller than HSI footprints.
                if ~isempty(mask_mastcam)
                    % nextplot_cur = obj.isvimage_MASTCAM_SelectMask.ax.NextPlot;
                    % obj.isvimage_MASTCAM_SelectMask.ax.NextPlot = 'add';
                    imobj = obj.update_isvimage_MASTCAM_SelectMask(mask_mastcam,[0.8,0.4,0]);
                    cursor_obj.UserData.MSTprojViewPlotObj.imObj_msldemc_UFOVmask_msldemresol = imobj;
                    % obj.isvimage_MASTCAM_SelectMask.ax.NextPlot = nextplot_cur;

                    lnObj_spc_mastcam = obj.plot_mastcam(obj.MSTview,mask_mastcam);
                    cursor_obj.UserData.MSTprojViewPlotObj.lnObj_spc_mastcam_msldemresol = lnObj_spc_mastcam;
                end
                if ~isempty(mask_msldemc_UFOVmask)
                    % nextplot_cur = obj.isvimage_NEE_SelectMask.ax.NextPlot;
                    % obj.isvimage_NEE_SelectMask.ax.NextPlot = 'add';
                    imobj = obj.update_isvimage_NEE_SelectMask(mask_msldemc_UFOVmask,[0.8,0.4,0]);
                    cursor_obj.UserData.MSTprojViewPlotObj.imObj_mastcam_msldemresol = imobj;
                    % obj.isvimage_NEE_SelectMask.ax.NextPlot = nextplot_cur;
                end
            end
            
        end
        
        function [out] = ISV_MASTCAM_WindowKeyPressFcn(obj,figobj,eventData)
            [out] = obj.ISV_MASTCAM.obj_HSIview.image_WindowKeyPressFcn_HSIview(figobj,eventData);
            % future implementation
        end
        
        function ISV_MASTCAM_BtnDwnFcn_update_mask_MSLDEMresol(obj,cursor_obj,x_mst,y_mst)
            [x_east,y_north,x_dem,y_dem,mask_mastcam,mask_msldemc_UFOVmask]...
                = obj.MSTproj.create_SelectMask_fromMASTCAM(x_mst,y_mst);
            % x_dem = obj.MSTproj.mastcam_msldemc_nn_UFOVmask(y_mst,x_mst,1);
            % y_dem = obj.MSTproj.mastcam_msldemc_nn_UFOVmask(y_mst,x_mst,2);

            % MASK for MSLDEM image
            % [mask_msldemc_UFOVmask] = obj.create_mask_msldem_wMSLDEMcoord(x_dem,y_dem,x_mst,y_mst);
            % show the mask on ImageStackView for MSLDEM in an independent
            % layer.
            obj.update_isvimage_NEE_SelectMask(mask_msldemc_UFOVmask);

            % With reference MSLDEM image coordinate (x_dem,y_dem), create
            % mask_mastcam (boolean image showing the corresponding pixels)
            % [mask_mastcam] = obj.create_mask_mastcam_wMSLDEMcoord(x_dem,y_dem,x_mst,y_mst);

            % show the mask on ImageStackView for MASTCAM in an independent
            % layer.
            obj.update_isvimage_MASTCAM_SelectMask(mask_mastcam);
            lnObj_spc_mastcam = obj.plot_mastcam(obj.MSTview,mask_mastcam);
            hsivplot_obj = HSIviewPlot();
            hsivplot_obj.cursor_obj = cursor_obj;
            hsivplot_obj.line_obj = lnObj_spc_mastcam;
            cursor_obj.UserData.HSIviewPlot_obj = [cursor_obj.UserData.HSIviewPlot_obj hsivplot_obj];
        end
        
        function [] = ISV_MASTCAM_BtnDwnFcn_HSI(obj,x_east,y_north,mask_msldemc_UFOVmask)
            for i=1:obj.objHSIview.nhsi
                hsielem = obj.objHSIview.hsiar(i);
                % get range of the selected hsi pixel.
                [mask_hsi,mask_msldemc_UFOVmask_hsiresol,...
                    mask_mastcam_hsiresol] = obj.get_hsi_mask(...
                    hsielem,x_east,y_north,mask_msldemc_UFOVmask);
                
                nextplot_cur = obj.isvimage_NEE_SelectMask.ax.NextPlot;
                obj.isvimage_NEE_SelectMask.ax.NextPlot = 'add';
                obj.update_isvimage_NEE_SelectMask(mask_msldemc_UFOVmask_hsiresol);
                obj.isvimage_NEE_SelectMask.ax.NextPlot = nextplot_cur;
                
                lnObj_spc_mastcam = obj.plot_mastcam(obj.MSTview,mask_mastcam);
                
                % hsii_i_NEEmask = false(hsiari.hdr.lines,hsiari.hdr.samples);
                % hsii_i_NEEmask(l_proj,s_proj) = true;
                
                obj.update_isvimage_MASTCAM_SelectMask(mask_mastcam_hsiresol);

                
            end
        end
        
        %%
        % callback functions utilities.      
        function [mask_hsi,mask_msldemc_UFOVmask_hsiresol,...
                mask_mastcam_hsiresol]...
                = get_hsi_mask(obj,hsielem,x_east,y_north,mask_msldemc_UFOVmask)
            mask_hsi = false(hsielem.imszy,hsielem.imszx);
            % nearest HSI image pixel.
            [xi,yi] = hsielem.get_xy_fromNE(x_east,y_north);
            mask_hsi(yi,xi) = true;
            [pos_row,pos_col] = find(mask_msldemc_UFOVmask);
            for i=1:length(pos_row)
                rowi = pos_row(i); coli = pos_col(i);
                [xei,yni] = obj.MSTproj.get_NE_from_msldemc_imUFOVxy(coli,rowi);
                % Get the HSI pixel nearest from one selected pixel in
                % msldemc_UFOVmask.
                [xi,yi] = hsielem.get_xy_fromNE(xei,yni);
                % mask_hsi_ar(i,1) = xi; mask_hsi_ar(i,2) = yi; 
                mask_hsi(yi,xi) = true;
            end
            mask_msldemc_UFOVmask_hsiresol = false(...
                obj.MSTproj.msldemc_imUFOVhdr.lines,...
                obj.MSTproj.msldemc_imUFOVhdr.samples);
            
            mask_mastcam_hsiresol = false(obj.MSTproj.MASTCAMdata.L_im,...
                obj.MSTproj.MASTCAMdata.S_im);
            [hsirow,hsicol] = find(mask_hsi);
            
            x1 = obj.MSTproj.msldemc_imUFOVhdr.x(1);
            dx = obj.MSTproj.MSLDEMdata.hdr.map_info.dx;
            y1 = obj.MSTproj.msldemc_imUFOVhdr.y(1);
            dy = obj.MSTproj.MSLDEMdata.hdr.map_info.dy;
            
            
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
                        if obj.MSTproj.msldemc_imUFOVmask(y_demi,x_demi)
                            idxes = obj.MSTproj.mapper_msldemc2mastcam_cell{obj.MSTproj.mapper_msldemc2mastcam_mat(y_demi,x_demi)};
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
        
        function [imobj] = update_isvimage_MASTCAM_SelectMask(obj,mask_mastcam,col)
            % show the mask on ImageStackView for MASTCAM in an
            % independent layer.
            col3 = reshape(col,[1,1,3]);
            if isempty(obj.isvimage_MASTCAM_SelectMask) || ~isvalid(obj.isvimage_MASTCAM_SelectMask)
                obj.isvimage_MASTCAM_SelectMask ...
                    = obj.MSTview.obj_HSIview.obj_ISV.add_layer(...
                    mask_mastcam.*col3,'AlphaData',double(mask_mastcam),'name','SelectMask');
                % obj.isvimage_MASTCAM_SelectMask.ax.NextPlot = 'replacechildren';
                obj.isvimage_MASTCAM_SelectMask.ismask = true;
                imobj = obj.isvimage_MASTCAM_SelectMask.imobj;
            else
                imobj = obj.isvimage_MASTCAM_SelectMask.add_children({...
                    mask_mastcam.*col3,'AlphaData',double(mask_mastcam)});
            end
        end
        
        function [imobj] = update_isvimage_NEE_SelectMask(obj,mask_msldemc_UFOVmask,col)
            col3 = reshape(col,[1,1,3]);
            if isempty(obj.isvimage_NEE_SelectMask) || ~isvalid(obj.isvimage_NEE_SelectMask)
                obj.isvimage_NEE_SelectMask = obj.ISV_MSLDEM.add_layer(...
                    obj.MSTproj.msldemc_imUFOVhdr.x([1,end]),...
                    obj.MSTproj.msldemc_imUFOVhdr.y([1,end]),...
                    mask_msldemc_UFOVmask.*col3,...
                    'AlphaData',double(mask_msldemc_UFOVmask),'name','SelectMask');
                obj.isvimage_NEE_SelectMask.ismask = true;
                imobj = obj.isvimage_NEE_SelectMask.imobj;
            else
                imobj = obj.isvimage_NEE_SelectMask.add_children(...
                    {obj.MSTproj.msldemc_imUFOVhdr.x([1,end]),...
                    obj.MSTproj.msldemc_imUFOVhdr.y([1,end]),...
                    mask_msldemc_UFOVmask.*col3,...
                    'AlphaData',double(mask_msldemc_UFOVmask)});
            end
        end
        
        function line_obj = plot_mastcam(obj,objMSTview,mask_mastcam,varargin)
            mask_mastcam_1nan = convertBoolTo1nan(mask_mastcam);
            spc = squeeze(nanmean(objMSTview.MSTMSI.img.*mask_mastcam_1nan,[1,2]));
            wv = objMSTview.MSTMSI.wavelength;
            bdxes = 1:objMSTview.MSTMSI.hdr.bands;
            line_obj = objMSTview.obj_HSIview.obj_SpecView.plot(...
                [wv,spc,'DisplayName',sprintf('MASTCAM'),varargin],...
                {'Band',bdxes});
            objMSTview.obj_HSIview.obj_SpecView.set_xlim();
            objMSTview.obj_HSIview.obj_SpecView.set_ylim();
        end
        
        function plot_NEE(obj,mask_NEE)
            
        end
        

    end
end

