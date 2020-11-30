classdef MASTCAMOrthoISV < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ISV
        MSTview
        CRISMview
    end
    
    methods
        function obj =  MASTCAMOrthoISV(mstview,varargin)
            %
            if (rem(length(varargin),2)==1)
                error('Optional parameters should always go by pairs');
            else
                for i=1:2:(length(varargin)-1)
                    switch upper(varargin{i})
                        % case {'SPECVIEW'}
                            % objSpecView = varargin{i+1};
                            % varargin_rmidx = [varargin_rmidx i i+1];
                        otherwise
                            error('Unrecognized option: %s',varargin{i});
                    end
                end
            end
            
            obj.MSTview = mstview;
            obj.init_ISV();
            
        end
        
        function init_ISV(obj)
            % get MSLDEM image
            % set margin pixels for the UFOVmask
            mrgn = 100;
            smpl_offset = max(1,obj.MSTview.MSTproj.msldemc_imUFOVhdr.sample_offset-mrgn);
            smpls = min(obj.MSTview.MSTproj.MSLDEMdata.hdr.samples-smpl_offset,...
                obj.MSTview.MSTproj.msldemc_imUFOVhdr.samples+2*mrgn);
            ln_offset = max(1,obj.MSTview.MSTproj.msldemc_imUFOVhdr.line_offset-mrgn);
            lns = min(obj.MSTview.MSTproj.MSLDEMdata.hdr.lines-ln_offset,...
                obj.MSTview.MSTproj.msldemc_imUFOVhdr.lines+2*mrgn);
            
            tic; msldemc_img = msldem_lazyenvireadRect(obj.MSTview.MSTproj.MSLDEMdata,...
                smpl_offset,ln_offset,smpls,lns,'precision','single'); toc;           
            
            obj.ISV.add_layer(...
                obj.MSTview.MSTproj.MSLDEMdata.hdr.x([smpl_offset+1,smpl_offset+smpls]),...
                obj.MSTview.MSTproj.MSLDEMdata.hdr.y([ln_offset+1,ln_offset+lns]),...
                msldemc_img,'name','MSLDEM');
            obj.ISV.add_layer(...
                obj.MSTview.MSTproj.msldemc_imUFOVhdr.x([1,end]),...
                obj.MSTview.MSTproj.msldemc_imUFOVhdr.y([1,end]),...
                obj.MSTview.MSTproj.msldemc_imUFOVmask,'name','UFOV mask',...
                'alphadata',obj.MSTview.MSTproj.msldemc_imUFOVmask>0);
            
            obj.ISV.Update_axim_aspectR();
            obj.ISV.Restore_ImageAxes2LimHome();

        end
        
        %%
        % =================================================================
        % ORTHO image button down callback function
        % =================================================================
        
        function [out] = MSTOrthoISV_BtnDwnFcn(obj,hObject,eventData)
            if ~isempty(obj.CRISMview) && isvalid(obj.CRISMview) && length(obj.CRISMview)==1
                [out] = obj.CRISMview.image_BtnDwnFcn_HSIview(hObject,eventData);
            else
                [out] = obj.ISV.image_BtnDwnFcn(hObject,eventData);
            end
            cursor_obj = out.cursor_obj;
            if isfield(cursor_obj.UserData,'MSTprojViewPlotObj')
            else
                cursor_obj.UserData.MSTprojViewPlotObj = MSTprojViewPlot();
            end
            
            % Obtain coordinate in the MSLDEM image pixel domain.
            
            x_east = out.cursor_obj.X; y_north = out.cursor_obj.Y;
            
            for i=1:length(obj.MSTview)
                mstview = obj.MSTview(i);
                % =============================================================
                % First get the selected point/area in the MASTCAM image at 
                % MSLDEM resolution.
                [x_mst,y_mst,x_dem,y_dem,mask_mastcam,mask_msldemc_UFOVmask]...
                    = mstview.MSTproj.get_SelectMask_fromNorthEast(x_east,y_north);
                
                if ~isempty(obj.CRISMview)
                    for hi=1:obj.CRISMview.nhsi
                        hsielem = obj.CRISMview.hsiar(hi);
                        % get range of the selected hsi pixel.
                        [mask_hsi,mask_msldemc_UFOVmask_hsiresol,...
                            mask_mastcam_hsiresol] = obj.get_hsi_mask(...
                            hsielem,x_east,y_north,mask_msldemc_UFOVmask);

                        [imobj] = obj.update_isvimage_NEE_SelectMask(mask_msldemc_UFOVmask_hsiresol,[0,0,0.8]);
                        cursor_obj.UserData.MSTprojViewPlotObj.imObj_msldemc_UFOVmask_hsiresol(hi) = imobj;
                        lnObj_spc_mastcam = obj.plot_mastcam(mask_mastcam_hsiresol);
                        cursor_obj.UserData.MSTprojViewPlotObj.lnObj_spc_mastcam_hsiresol(hi) = lnObj_spc_mastcam;

                        [imobj] = mstview.update_isvimage_MASTCAM_SelectMask(mask_mastcam_hsiresol,[0,0,0.8]);
                        cursor_obj.UserData.MSTprojViewPlotObj.imObj_mastcam_hsiresol(hi) = imobj;

                    end
                end
                
                % Show mask at MSLDEM resolution.
                % This is added after HSI image mask as we expect this would be
                % smaller than HSI footprints.
                if ~isempty(mask_mastcam)
                    % nextplot_cur = obj.isvimage_MASTCAM_SelectMask.ax.NextPlot;
                    % obj.isvimage_MASTCAM_SelectMask.ax.NextPlot = 'add';
                    imobj = mstview.update_isvimage_MASTCAM_SelectMask(mask_mastcam,[0.8,0.4,0]);
                    cursor_obj.UserData.MSTprojViewPlotObj.imObj_msldemc_UFOVmask_msldemresol = imobj;
                    % obj.isvimage_MASTCAM_SelectMask.ax.NextPlot = nextplot_cur;

                    lnObj_spc_mastcam = mstview.plot_mastcam(mask_mastcam);
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
        
        function [out] = MSTOrthoISV_WindowKeyPressFcn(obj,figobj,eventData)
            [out] = obj.ISV_ORTHO.ISVWindowKeyPressFcn(figobj,eventData);
            % future implementation
        end
        
        % function [] = MSTOrthoISV_BtnDwnFcn_HSI(obj,cursor_obj,x_east,y_north,mask_msldemc_UFOVmask)
            
        % end
        
        %%
        
        function [imobj] = update_isvimage_NEE_SelectMask(obj,mask_msldemc_UFOVmask,col)
            col3 = reshape(col,[1,1,3]);
            if isempty(obj.isvimage_NEE_SelectMask) || ~isvalid(obj.isvimage_NEE_SelectMask)
                obj.isvimage_NEE_SelectMask = obj.ISV_ORTHO.add_layer(...
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
        
        function line_obj = plot_mastcam(obj,mask_mastcam,varargin)
            mask_mastcam_1nan = convertBoolTo1nan(mask_mastcam);
            for i=1:obj.obj_HSIview.nhsi
                mstmsi = obj.obj_HSIview.hsiar(i).hsi;
                spc = squeeze(nanmean(mstmsi.img.*mask_mastcam_1nan,[1,2]));
                wv = mstmsi.wavelength;
                bdxes = 1:mstmsi.hdr.bands;
                line_obj = obj.obj_HSIview.obj_SpecView.plot(...
                    [wv,spc,'DisplayName',...
                    sprintf('%s AVE',obj.obj_HSIview.hsiar(i).name),varargin],...
                    {'Band',bdxes});
            end
            obj.obj_HSIview.obj_SpecView.set_xlim();
            obj.obj_HSIview.obj_SpecView.set_ylim();
        end
        
        %%
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
        
        
    end
end

