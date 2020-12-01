classdef MASTCAMOrthoISV < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ISV
        MSTview
        CRISMview
        isvimage_NEE_SelectMask
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
            obj.ISV = ImageStackView({},varargin{:},'Ydir','normal',...
                    'XY_COORDINATE_SYSTEM','NorthEast',...
                    'IMAGE_CURSOR_FCN',@obj.MSTOrthoISV_BtnDwnFcn,...
                    'IMAGE_WINDOWKEYPRESS_FCN',@obj.MSTOrthoISV_WindowKeyPressFcn);    
        end
        
        
        
        function [crismview] = add_crism(obj,varargin)
            crismview = HSIview(varargin{:},'ISV',obj.ISV,...
                'VARARGIN_ImageStackView',{'IMAGE_CURSOR_FCN',@obj.MSTOrthoISV_BtnDwnFcn,...
                'IMAGE_WINDOWKEYPRESS_FCN',@obj.MSTOrthoISV_WindowKeyPressFcn});
            obj.CRISMview = crismview;
        end
        
        function add_mstview(obj,mstview)
            obj.MSTview = [obj.MSTview mstview];
            mstview.MSTOrthoISV = obj;
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
            
            % first create individual masks and combined masks.
            [mask_mastcam_cellar,mask_msldemc_UFOVmask_cellar,...
                msldemc_mask_all,msldemc_mask_all_hdr]...
                = mstproj_combine_mask([obj.MSTview.MSTproj],x_east,y_north);
            
            % if ~isempty(obj.CRISMview)
            %     obj.MSTOrthoISV_BtnDwnFcn_HSI(obj,cursor_obj,x_east,y_north,msldemc_mask_all);
            % end
            
            % finally update masks 
            for i=1:length(obj.MSTview)
                mstview = obj.MSTview(i);
                mstproj = mstview.MSTproj;
                si1 = mstproj.msldemc_imUFOVhdr.sample_offset-msldemc_mask_all_hdr.sample_offset+1;
                li1 = mstproj.msldemc_imUFOVhdr.line_offset-msldemc_mask_all_hdr.line_offset+1;
                siend = si1+mstproj.msldemc_imUFOVhdr.samples-1;
                liend = li1+mstproj.msldemc_imUFOVhdr.lines-1;
                [mask_mastcam_cmb,mask_msldemc_UFOVmask_cmb_bp]...
                    = mstproj.create_SelectMask_from_msldemc_mask(msldemc_mask_all(li1:liend,si1:siend));
                % =============================================================
                % First get the selected point/area in the MASTCAM image at 
                % MSLDEM resolution.
                [x_mst,y_mst,x_dem,y_dem,mask_mastcam,mask_msldemc_UFOVmask]...
                    = mstview.MSTproj.get_SelectMask_fromNorthEast(x_east,y_north);
                
                if ~isempty(obj.CRISMview)
                    obj.MSTOrthoISV_BtnDwnFcn_HSI(cursor_obj,mstview,mstproj,...
                        x_east,y_north,mask_msldemc_UFOVmask_cmb_bp);
                end
                
                
                % Show mask at MSLDEM resolution.
                % This is added after HSI image mask as we expect this would be
                % smaller than HSI footprints.
                if ~isempty(mask_mastcam_cmb)
                    % nextplot_cur = obj.isvimage_MASTCAM_SelectMask.ax.NextPlot;
                    % obj.isvimage_MASTCAM_SelectMask.ax.NextPlot = 'add';
                    imobj = mstview.update_isvimage_MASTCAM_SelectMask(mask_mastcam_cmb,[0.8,0.4,0]);
                    cursor_obj.UserData.MSTprojViewPlotObj.imObj_msldemc_UFOVmask_msldemresol = [cursor_obj.UserData.MSTprojViewPlotObj.imObj_msldemc_UFOVmask_msldemresol imobj];
                    % obj.isvimage_MASTCAM_SelectMask.ax.NextPlot = nextplot_cur;

                    lnObj_spc_mastcam = mstview.plot_mastcam(mask_mastcam_cmb);
                    cursor_obj.UserData.MSTprojViewPlotObj.lnObj_spc_mastcam_msldemresol = [cursor_obj.UserData.MSTprojViewPlotObj.lnObj_spc_mastcam_msldemresol lnObj_spc_mastcam];
                end
                if ~isempty(mask_msldemc_UFOVmask_cmb_bp)
                    % nextplot_cur = obj.isvimage_NEE_SelectMask.ax.NextPlot;
                    % obj.isvimage_NEE_SelectMask.ax.NextPlot = 'add';
                    imobj = obj.update_isvimage_NEE_SelectMask(mstproj,mask_msldemc_UFOVmask_cmb_bp,[0.8,0.4,0]);
                    cursor_obj.UserData.MSTprojViewPlotObj.imObj_mastcam_msldemresol = [cursor_obj.UserData.MSTprojViewPlotObj.imObj_mastcam_msldemresol imobj];
                    % obj.isvimage_NEE_SelectMask.ax.NextPlot = nextplot_cur;
                end
            
            end
            
        
        end
        
        function [out] = MSTOrthoISV_WindowKeyPressFcn(obj,figobj,eventData)
            [out] = obj.ISV.ISVWindowKeyPressFcn(figobj,eventData);
            % future implementation
        end
        
        function [] = MSTOrthoISV_BtnDwnFcn_HSI(obj,cursor_obj,mstview,mstproj,x_east,y_north,mask_msldemc_UFOVmask)
            for hi=1:obj.CRISMview.nhsi
                hsielem = obj.CRISMview.hsiar(hi);
                % get range of the selected hsi pixel.
                [mask_hsi,mask_msldemc_UFOVmask_hsiresol,...
                    mask_mastcam_hsiresol] = mstproj_get_hsi_mask(mstproj,...
                    hsielem,x_east,y_north,mask_msldemc_UFOVmask);

                [imobj] = obj.update_isvimage_NEE_SelectMask(mstproj,mask_msldemc_UFOVmask_hsiresol,[0,0,0.8]);
                cursor_obj.UserData.MSTprojViewPlotObj.imObj_msldemc_UFOVmask_hsiresol = [cursor_obj.UserData.MSTprojViewPlotObj.imObj_msldemc_UFOVmask_hsiresol imobj];
                lnObj_spc_mastcam = mstview.plot_mastcam(mask_mastcam_hsiresol);
                cursor_obj.UserData.MSTprojViewPlotObj.lnObj_spc_mastcam_hsiresol = [cursor_obj.UserData.MSTprojViewPlotObj.lnObj_spc_mastcam_hsiresol lnObj_spc_mastcam];

                [imobj] = mstview.update_isvimage_MASTCAM_SelectMask(mask_mastcam_hsiresol,[0,0,0.8]);
                cursor_obj.UserData.MSTprojViewPlotObj.imObj_mastcam_hsiresol = [cursor_obj.UserData.MSTprojViewPlotObj.imObj_mastcam_hsiresol imobj];

            end
        end
        
        %%
        function MSTOrthoISV_mastcam(obj,cursor_obj,msldemc_mask_ref,mstproj_ref,x_east,y_north)
            
            [msldemc_mask_all,msldemc_mask_all_hdr]...
                = mstproj_get_msldem_mask_from_msldem_mask([obj.MSTview.MSTproj],msldemc_mask_ref,mstproj_ref);
            
            for i=1:length(obj.MSTview)
                mstview = obj.MSTview(i);
                mstproj = mstview.MSTproj;
                si1 = mstproj.msldemc_imUFOVhdr.sample_offset-msldemc_mask_all_hdr.sample_offset+1;
                li1 = mstproj.msldemc_imUFOVhdr.line_offset-msldemc_mask_all_hdr.line_offset+1;
                siend = si1+mstproj.msldemc_imUFOVhdr.samples-1;
                liend = li1+mstproj.msldemc_imUFOVhdr.lines-1;
                [mask_mastcam_cmb,mask_msldemc_UFOVmask_cmb_bp]...
                    = mstproj.create_SelectMask_from_msldemc_mask(msldemc_mask_all(li1:liend,si1:siend));
                % =============================================================
                % First get the selected point/area in the MASTCAM image at 
                % MSLDEM resolution.
                
                if ~isempty(obj.CRISMview)
                    obj.MSTOrthoISV_BtnDwnFcn_HSI(cursor_obj,mstview,mstproj,...
                        x_east,y_north,mask_msldemc_UFOVmask_cmb_bp);
                end
                
                
                % Show mask at MSLDEM resolution.
                % This is added after HSI image mask as we expect this would be
                % smaller than HSI footprints.
                if ~isempty(mask_mastcam_cmb)
                    % nextplot_cur = obj.isvimage_MASTCAM_SelectMask.ax.NextPlot;
                    % obj.isvimage_MASTCAM_SelectMask.ax.NextPlot = 'add';
                    imobj = mstview.update_isvimage_MASTCAM_SelectMask(mask_mastcam_cmb,[0.8,0.4,0]);
                    cursor_obj.UserData.MSTprojViewPlotObj.imObj_msldemc_UFOVmask_msldemresol = [cursor_obj.UserData.MSTprojViewPlotObj.imObj_msldemc_UFOVmask_msldemresol imobj];
                    % obj.isvimage_MASTCAM_SelectMask.ax.NextPlot = nextplot_cur;

                    lnObj_spc_mastcam = mstview.plot_mastcam(mask_mastcam_cmb);
                    cursor_obj.UserData.MSTprojViewPlotObj.lnObj_spc_mastcam_msldemresol = [cursor_obj.UserData.MSTprojViewPlotObj.lnObj_spc_mastcam_msldemresol lnObj_spc_mastcam];
                end
                if ~isempty(mask_msldemc_UFOVmask_cmb_bp)
                    % nextplot_cur = obj.isvimage_NEE_SelectMask.ax.NextPlot;
                    % obj.isvimage_NEE_SelectMask.ax.NextPlot = 'add';
                    imobj = obj.update_isvimage_NEE_SelectMask(mstproj,mask_msldemc_UFOVmask_cmb_bp,[0.8,0.4,0]);
                    cursor_obj.UserData.MSTprojViewPlotObj.imObj_mastcam_msldemresol = [cursor_obj.UserData.MSTprojViewPlotObj.imObj_mastcam_msldemresol imobj];
                    % obj.isvimage_NEE_SelectMask.ax.NextPlot = nextplot_cur;
                end
            
            end
            
        end
        
        %%
        
        function [imobj] = update_isvimage_NEE_SelectMask(obj,mstproj,mask_msldemc_UFOVmask,col)
            col3 = reshape(col,[1,1,3]);
            if isempty(obj.isvimage_NEE_SelectMask) || ~isvalid(obj.isvimage_NEE_SelectMask)
                obj.isvimage_NEE_SelectMask = obj.ISV.add_layer(...
                    mstproj.msldemc_imUFOVhdr.x([1,end]),...
                    mstproj.msldemc_imUFOVhdr.y([1,end]),...
                    mask_msldemc_UFOVmask.*col3,...
                    'AlphaData',double(mask_msldemc_UFOVmask),'name','SelectMask');
                obj.isvimage_NEE_SelectMask.ismask = true;
                imobj = obj.isvimage_NEE_SelectMask.imobj;
            else
                imobj = obj.isvimage_NEE_SelectMask.add_children(...
                    {mstproj.msldemc_imUFOVhdr.x([1,end]),...
                    mstproj.msldemc_imUFOVhdr.y([1,end]),...
                    mask_msldemc_UFOVmask.*col3,...
                    'AlphaData',double(mask_msldemc_UFOVmask)});
            end
        end
        
    end
end

