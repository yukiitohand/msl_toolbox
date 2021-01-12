classdef MASTCAMMSIview < handle
    %MASTCAMview
    %   Spectral Viewer for MASTCAM images
    properties
        obj_HSIview
        MASTCAMdataseq_eye
        MSTMSI
        MSTproj
        MSTOrthoISV
        % isvimage_NEE_SelectMask
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
            
        function add_projection(obj,MSTproj,varargin)
            if (rem(length(varargin),2)==1)
                error('Optional parameters should always go by pairs');
            else
                for i=1:2:(length(varargin)-1)
                    switch upper(varargin{i})
                        case {'MASTCAMORTHOVIEW'}
                            obj.MSTOrthoISV = varargin{i+1};
                            % varargin_rmidx = [varargin_rmidx i i+1];
                        otherwise
                            error('Unrecognized option: %s',varargin{i});
                    end
                end
            end
            obj.obj_HSIview.obj_ISV.add_layer(MSTproj.mastcam_range,'name','Range');
            obj.obj_HSIview.obj_ISV.add_layer(MSTproj.mastcam_emi,'name','Emission');
            obj.MSTproj = MSTproj;
            if isempty(obj.MSTOrthoISV)
                obj.MSTOrthoISV = MASTCAMOrthoISV(obj);
            else
                obj.MSTOrthoISV.add_mstview(obj);
            end
            obj.init_MSTOrthoISV();
            
        end
        
        function init_MSTOrthoISV(obj)
            % get MSLDEM image
            % set margin pixels for the UFOVmask
            [msldemc_img,msldemc_img_hdr] = obj.MSTproj.MSLDEMread();
            
            obj.MSTOrthoISV.ISV.add_layer(...
                obj.MSTproj.MSLDEMdata.hdr.x([msldemc_img_hdr.sample_offset+1,msldemc_img_hdr.sample_offset+msldemc_img_hdr.samples]),...
                obj.MSTproj.MSLDEMdata.hdr.y([msldemc_img_hdr.line_offset+1,msldemc_img_hdr.line_offset+msldemc_img_hdr.lines]),...
                msldemc_img,'name','MSLDEM');
            obj.MSTOrthoISV.ISV.add_layer(...
                obj.MSTproj.msldemc_imUFOVhdr.x([1,end]),...
                obj.MSTproj.msldemc_imUFOVhdr.y([1,end]),...
                obj.MSTproj.msldemc_imUFOVmask,'name','UFOV mask',...
                'alphadata',obj.MSTproj.msldemc_imUFOVmask>0);
            
            obj.MSTOrthoISV.ISV.Update_axim_aspectR();
            obj.MSTOrthoISV.ISV.Restore_ImageAxes2LimHome();

        end
        
        %%
        % =================================================================
        % ORTHO image button down callback function
        % =================================================================
        
        
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
            
            
            if ~isempty(obj.MSTOrthoISV) && isvalid(obj.MSTOrthoISV)
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

                    if ~isempty(obj.MSTOrthoISV)
                        obj.MSTOrthoISV.MSTOrthoISV_mastcam(cursor_obj,mask_msldemc_UFOVmask,obj.MSTproj,x_east,y_north);
                    end

                    % Show mask at MSLDEM resolution.
                    % This is added after HSI image mask as we expect this would be
                    % smaller than HSI footprints.
                    if ~isempty(mask_mastcam)
                        % nextplot_cur = obj.isvimage_MASTCAM_SelectMask.ax.NextPlot;
                        % obj.isvimage_MASTCAM_SelectMask.ax.NextPlot = 'add';
                        imobj = obj.update_isvimage_MASTCAM_SelectMask(mask_mastcam,[0.8,0.4,0]);
                        cursor_obj.UserData.MSTprojViewPlotObj.imObj_msldemc_UFOVmask_msldemresol = [cursor_obj.UserData.MSTprojViewPlotObj.imObj_msldemc_UFOVmask_msldemresol imobj];
                        % obj.isvimage_MASTCAM_SelectMask.ax.NextPlot = nextplot_cur;

                        lnObj_spc_mastcam = obj.plot_mastcam(mask_mastcam);
                        cursor_obj.UserData.MSTprojViewPlotObj.lnObj_spc_mastcam_msldemresol = [cursor_obj.UserData.MSTprojViewPlotObj.lnObj_spc_mastcam_msldemresol lnObj_spc_mastcam];
                    end
                    if ~isempty(mask_msldemc_UFOVmask)
                        % nextplot_cur = obj.isvimage_NEE_SelectMask.ax.NextPlot;
                        % obj.isvimage_NEE_SelectMask.ax.NextPlot = 'add';
                        imobj = obj.MSTOrthoISV.update_isvimage_NEE_SelectMask(obj.MSTproj,mask_msldemc_UFOVmask,[0.8,0.4,0]);
                        cursor_obj.UserData.MSTprojViewPlotObj.imObj_mastcam_msldemresol = [cursor_obj.UserData.MSTprojViewPlotObj.imObj_mastcam_msldemresol imobj];
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
        
        
        
    end
end