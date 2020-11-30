classdef MASTCAMMSIview < handle
    %MASTCAMview
    %   Spectral Viewer for MASTCAM images
    properties
        obj_HSIview
        MASTCAMdataseq_eye
        MSTMSI
        MSTproj
        ISV_ORTHO
        isvimage_NEE_SelectMask
        isvimage_MASTCAM_SelectMask
        CRISMview
        MSTMSIview
    end
    
    methods
        function obj = MASTCAMMSIview(mstmsi_input,varargin)
            objSpecView = [];
            % varargin_rmidx = [];
            if (rem(length(varargin),2)==1)
                error('Optional parameters should always go by pairs');
            else
                for i=1:2:(length(varargin)-1)
                    switch upper(varargin{i})
                        case {'SPECVIEW'}
                            objSpecView = varargin{i+1};
                            % varargin_rmidx = [varargin_rmidx i i+1];
                        otherwise
                            error('Unrecognized option: %s',varargin{i});
                    end
                end
            end
            % varargin_HSI = varargin(setdiff(1:length(varargin),varargin_rmidx));
            %Constructor
            % obj.MSTMSI = mstmsi;
            
            
            % for i=1:length(hsiar_input)
            [rgbList] = obj.process_mstmsi_input(mstmsi_input);
            if ~iscell(rgbList),rgbList = {rgbList}; end
            obj.obj_HSIview = HSIview(rgbList,...
            mstmsi_input,...
            ...'SPC_XLIM',[300 1200],...
            'varargin_ImageStackView',{'Ydir','reverse','XY_COORDINATE_SYSTEM','IMAGEPIXELS'},...
            'SpecView',objSpecView);
        
            obj.obj_HSIview.obj_ISV.custom_image_cursor_fcn = @obj.ISV_MASTCAM_BtnDwnFcn;
            obj.obj_HSIview.obj_ISV.custom_windowkeypress_fcn = @obj.ISV_MASTCAM_WindowKeyPressFcn;
            
        end
        
        function [rgbList] = process_mstmsi_input(obj,mstmsi_input)

            if isempty(mstmsi_input)
                rgbList = [];
                % skip if input hsiar is empty.
            elseif isa(mstmsi_input,'MASTCAMMSI')
                rgbList = mstmsi_input.get_rgb();
            elseif iscell(mstmsi_input)
                if length(mstmsi_input)==1
                    if isa(mstmsi_input{1},'MASTCAMMSI')
                        rgbList = mstmsi_input{1}.get_rgb();
                    elseif iscell(mstmsi_input{1}) && isa(mstmsi_input{1}{1},'MASTCAMMSI')
                        rgbList = mstmsi_input{1}{1}.get_rgb();
                    else
                        error('Input is not proper.');
                    end
                elseif length(mstmsi_input)>1
                    if isa(mstmsi_input{1},'MASTCAMMSI') && ~(isa(mstmsi_input{2},'MASTCAMMSI') || iscell(mstmsi_input{2}))
                        rgbList = mstmsi_input{1}.get_rgb();
                    else
                        nelem = length(mstmsi_input);
                        rgbList = cell(1,nelem);
                        for i=1:nelem
                            if iscell(mstmsi_input{i})
                                rgbList{i} = mstmsi_input{i}{1}.get_rgb();
                            elseif isa(mstmsi_input{i},'MASTCAMMSI')
                                rgbList{i} = mstmsi_input{i}.get_rgb();
                            else
                                error('Input mstmsi is not proper.');
                            end
                        end
                    end
                end
            end
        end
            
        function add_projection(obj,MSTproj)
            obj.obj_HSIview.obj_ISV.add_layer(MSTproj.mastcam_range,'name','Range');
            obj.MSTproj = MSTproj;
            if isempty(obj.ISV_ORTHO)
                obj.ISV_ORTHO = ImageStackView({},'Ydir','normal',...
                    'XY_COORDINATE_SYSTEM','NorthEast',...
                    'IMAGE_CURSOR_FCN',@obj.ISV_ORTHO_BtnDwnFcn,...
                    'IMAGE_WINDOWKEYPRESS_FCN',@obj.ISV_ORTHO_WindowKeyPressFcn);
            end
            obj.init_ISV_ORTHO();
            
            % if isempty(obj.objHSIview)
            %     
            % else
            %     obj.ISV_MSLDEM = obj.objHSIview.obj_ISV;
            %     obj.ISV_MSLDEM.custom_image_cursor_fcn = @obj.ISV_MSLDEM_BtnDwnFcn;
            %     obj.ISV_MSLDEM.custom_windowkeypress_fcn = @obj.ISV_MSLDEM_WindowKeyPressFcn;
            %     obj.init_ISV_MSLDEM();
            % end
            
        end
        
        function init_ISV_ORTHO(obj)
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
            
            obj.ISV_ORTHO.add_layer(...
                obj.MSTproj.MSLDEMdata.hdr.x([smpl_offset+1,smpl_offset+smpls]),...
                obj.MSTproj.MSLDEMdata.hdr.y([ln_offset+1,ln_offset+lns]),...
                msldemc_img,'name','MSLDEM');
            obj.ISV_ORTHO.add_layer(...
                obj.MSTproj.msldemc_imUFOVhdr.x([1,end]),...
                obj.MSTproj.msldemc_imUFOVhdr.y([1,end]),...
                obj.MSTproj.msldemc_imUFOVmask,'name','UFOV mask',...
                'alphadata',obj.MSTproj.msldemc_imUFOVmask>0);
            
            obj.ISV_ORTHO.Update_axim_aspectR();
            obj.ISV_ORTHO.Restore_ImageAxes2LimHome();

        end
        
        function [crismview] = add_crism(obj,varargin)
            crismview = HSIview(varargin{:},'ISV',obj.ISV_ORTHO,...
                'SpecView',obj.obj_HSIview.obj_SpecView,...
                'VARARGIN_ImageStackView',{'IMAGE_CURSOR_FCN',@obj.ISV_ORTHO_BtnDwnFcn,...
                'IMAGE_WINDOWKEYPRESS_FCN',@obj.ISV_ORTHO_WindowKeyPressFcn});
            obj.CRISMview = crismview;
        end
        
        function [mstview] = add_mastcammsi(obj,varargin)
            mstview = MASTCAMMSIview(varargin{:},'ISV',obj.ISV_ORTHO,...
                'SpecView',obj.obj_HSIview.obj_SpecView,...
                'VARARGIN_ImageStackView',{'IMAGE_CURSOR_FCN',@obj.ISV_ORTHO_BtnDwnFcn,...
                'IMAGE_WINDOWKEYPRESS_FCN',@obj.ISV_ORTHO_WindowKeyPressFcn});
            obj.MSTMSIview = mstview;
        end
        
        %%
        % =================================================================
        % ORTHO image button down callback function
        % =================================================================
        
        function [out] = ISV_ORTHO_BtnDwnFcn(obj,hObject,eventData)
            if ~isempty(obj.CRISMview) && isvalid(obj.CRISMview)
                [out] = obj.CRISMview.image_BtnDwnFcn_HSIview(hObject,eventData);
            else
                [out] = obj.ISV_ORTHO.image_BtnDwnFcn(hObject,eventData);
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
            if ~isempty(obj.CRISMview)
                obj.ISV_ORTHO_BtnDwnFcn_HSI(cursor_obj,x_east,y_north,...
                    mask_msldemc_UFOVmask);
            end
            
            % Show mask at MSLDEM resolution.
            % This is added after HSI image mask as we expect this would be
            % smaller than HSI footprints.
            if ~isempty(mask_mastcam)
                % nextplot_cur = obj.isvimage_MASTCAM_SelectMask.ax.NextPlot;
                % obj.isvimage_MASTCAM_SelectMask.ax.NextPlot = 'add';
                imobj = obj.update_isvimage_MASTCAM_SelectMask(mask_mastcam,[0.8,0.4,0]);
                cursor_obj.UserData.MSTprojViewPlotObj.imObj_msldemc_UFOVmask_msldemresol = imobj;
                % obj.isvimage_MASTCAM_SelectMask.ax.NextPlot = nextplot_cur;
                
                lnObj_spc_mastcam = obj.plot_mastcam(mask_mastcam);
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
        
        function [out] = ISV_ORTHO_WindowKeyPressFcn(obj,figobj,eventData)
            [out] = obj.ISV_ORTHO.ISVWindowKeyPressFcn(figobj,eventData);
            % future implementation
        end
        
        function [] = ISV_ORTHO_BtnDwnFcn_HSI(obj,cursor_obj,x_east,y_north,mask_msldemc_UFOVmask)
            for i=1:obj.CRISMview.nhsi
                hsielem = obj.CRISMview.hsiar(i);
                % get range of the selected hsi pixel.
                [mask_hsi,mask_msldemc_UFOVmask_hsiresol,...
                    mask_mastcam_hsiresol] = obj.get_hsi_mask(...
                    hsielem,x_east,y_north,mask_msldemc_UFOVmask);

                [imobj] = obj.update_isvimage_NEE_SelectMask(mask_msldemc_UFOVmask_hsiresol,[0,0,0.8]);
                cursor_obj.UserData.MSTprojViewPlotObj.imObj_msldemc_UFOVmask_hsiresol(i) = imobj;
                lnObj_spc_mastcam = obj.plot_mastcam(mask_mastcam_hsiresol);
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
            [out] = obj.obj_HSIview.image_BtnDwnFcn_HSIview(hObject,eventData);
            cursor_obj = out.cursor_obj;
            if isfield(cursor_obj.UserData,'MSTprojViewPlotObj')
            else
                cursor_obj.UserData.MSTprojViewPlotObj = MSTprojViewPlot();
            end
            
            
            if ~isempty(obj.ISV_ORTHO) && isvalid(obj.ISV_ORTHO)
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
                    obj.ISV_ORTHO_BtnDwnFcn_HSI(cursor_obj,x_east,y_north,mask_msldemc_UFOVmask);

                    % Show mask at MSLDEM resolution.
                    % This is added after HSI image mask as we expect this would be
                    % smaller than HSI footprints.
                    if ~isempty(mask_mastcam)
                        % nextplot_cur = obj.isvimage_MASTCAM_SelectMask.ax.NextPlot;
                        % obj.isvimage_MASTCAM_SelectMask.ax.NextPlot = 'add';
                        imobj = obj.update_isvimage_MASTCAM_SelectMask(mask_mastcam,[0.8,0.4,0]);
                        cursor_obj.UserData.MSTprojViewPlotObj.imObj_msldemc_UFOVmask_msldemresol = imobj;
                        % obj.isvimage_MASTCAM_SelectMask.ax.NextPlot = nextplot_cur;

                        lnObj_spc_mastcam = obj.plot_mastcam(mask_mastcam);
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
            
        end
        
        function [out] = ISV_MASTCAM_WindowKeyPressFcn(obj,figobj,eventData)
            [out] = obj.ISV_MASTCAM.obj_HSIview.image_WindowKeyPressFcn_HSIview(figobj,eventData);
            % future implementation
        end
        
        %%
        function [imobj] = update_isvimage_MASTCAM_SelectMask(obj,mask_mastcam,col)
            % show the mask on ImageStackView for MASTCAM in an
            % independent layer.
            col3 = reshape(col,[1,1,3]);
            if isempty(obj.isvimage_MASTCAM_SelectMask) || ~isvalid(obj.isvimage_MASTCAM_SelectMask)
                obj.isvimage_MASTCAM_SelectMask ...
                    = obj.obj_HSIview.obj_ISV.add_layer(...
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
                                idxes_jj = obj.MSTproj.mapper_mastcam2msldemc{idxes(2,ii),idxes(1,ii)};
                                for jj=1:size(idxes_jj,2)
                                    mask_msldemc_UFOVmask_hsiresol(idxes_jj(2,jj),idxes_jj(1,jj)) = true;
                                end
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