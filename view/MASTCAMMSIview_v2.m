classdef MASTCAMMSIview_v2 < handle
    %MASTCAMview
    %   Spectral Viewer for MASTCAM images
    properties
        objENVIRasterMultview
        % MASTCAMdataseq_eye
        MSTMSI
        MSTproj
        objISV_proj
        % isvimage_NEE_SelectMask
        objISVImage_ISV_proj_SelectMask
        objISVImage_MASTCAM_SelectMask
    end
    
    methods
        function obj = MASTCAMMSIview_v2(mstmsi_input,varargin)
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
            obj.objENVIRasterMultview = ENVIRasterMultview(rgbList,...
            mstmsi_input,...
            ...'SPC_XLIM',[300 1200],...
            'varargin_ImageStackView',{'Ydir','reverse','XY_COORDINATE_SYSTEM','IMAGEPIXELS'},...
            'SpecView',objSpecView);
            obj.objENVIRasterMultview.obj_ISV.custom_image_cursor_fcn = @obj.ISV_MASTCAM_BtnDwnFcn;
            obj.objENVIRasterMultview.obj_ISV.custom_windowkeypress_fcn = @obj.ISV_MASTCAM_WindowKeyPressFcn;
            
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
            varargin_ISV = {};
            if (rem(length(varargin),2)==1)
                error('Optional parameters should always go by pairs');
            else
                for i=1:2:(length(varargin)-1)
                    switch upper(varargin{i})
                        case {'ISV_PROJ'}
                            obj.objISV_proj = varargin{i+1};
                        case {'ISV_VARARGIN'}
                            varargin_ISV = varargin{i+1};
                        otherwise
                            error('Unrecognized option: %s',varargin{i});
                    end
                end
            end
            % Add projection components
            mst_range_img = MSTproj.mastcam_proj.range.readimg('precision','raw');
            obj.objENVIRasterMultview.obj_ISV.add_layer(mst_range_img,'name','Range');
            mst_emiang_img = MSTproj.mastcam_proj.emiang.readimg('precision','raw');
            obj.objENVIRasterMultview.obj_ISV.add_layer(mst_emiang_img,'name','Emission');
            obj.MSTproj = MSTproj;
            
            if isempty(obj.objISV_proj)
                obj.objISV_proj = ImageStackView([], ...
                    'XY_COORDINATE_SYSTEM','LATLON', ...
                    'IMAGE_CURSOR_FCN',@obj.ISV_proj_BtnDwnFcn, ...
                    'IMAGE_WINDOWKEYPRESS_FCN',@obj.ISV_proj_WindowKeyPressFcn, ...
                    varargin_ISV{:});
            end
            
            % Show UFOV mask.
            msldemc_ufovmask_img = MSTproj.msldemc_imUFOV.mask.readimg('precision','raw');
            srng = MSTproj.msldemc_imUFOV.mask.get_lon_ctrrange();
            lrng = MSTproj.msldemc_imUFOV.mask.get_lat_ctrrange();
            isvimgObj_mask = obj.objISV_proj.add_layer(srng,lrng,msldemc_ufovmask_img,'name','ufov mask');
            isvimgObj_mask.imobj.AlphaData = double(msldemc_ufovmask_img>0)*0.5;
            obj.objISV_proj.Update_axim_aspectR();
            obj.objISV_proj.Restore_ImageAxes2LimHome();
        end
        
        %%
        % =================================================================
        % MASTCAM window button down callback function
        % =================================================================
        function create_ax_ISV_proj_SelectMask(obj)
            obj.objISVImage_ISV_proj_SelectMask = obj.objISV_proj.add_layer([],'name','Select Mask');
            obj.objISVImage_ISV_proj_SelectMask.ax.NextPlot = 'add';
        end
        
        function [out] = ISV_MASTCAM_BtnDwnFcn(obj,objAxes,eventData)
            [out] = obj.objENVIRasterMultview.image_BtnDwnFcn_HSIview(objAxes,eventData);
            cursor_obj = out.cursor_obj;
            if isfield(cursor_obj.UserData,'MSTprojViewPlotObj')
            else
                cursor_obj.UserData.MSTprojViewPlotObj = HSIviewPlot();
            end
            obj.ISV_MASTCAM_update_cursor(cursor_obj);
        end
        
        function [out] = ISV_MASTCAM_WindowKeyPressFcn(obj,figobj,eventData)
            [out] = obj.objENVIRasterMultview.image_WindowKeyPressFcn_HSIview(figobj,eventData);
            cursor_obj = out.cursor_obj;
            switch eventData.Key
                case {'rightarrow','leftarrow','uparrow','downarrow'}
                    obj.ISV_MASTCAM_update_cursor(cursor_obj);
            end
            
        end
        
        function [] = ISV_MASTCAM_update_cursor(obj,cursor_obj)
            
            if ~isempty(obj.objISV_proj) && isvalid(obj.objISV_proj)
                x_mst = cursor_obj.X; y_mst = cursor_obj.Y;
                x_mst_int = round(x_mst); y_mst_int = round(y_mst);
                % if cursor is pointing the sky, then do not any further
                % processing
                [x_demc,y_demc] = obj.MSTproj.mapper.get_msldemc_index(x_mst_int,y_mst_int);
                if isempty(x_demc)
                    if ~isempty(cursor_obj.UserData.MSTprojViewPlotObj.im_obj)
                        imObj = cursor_obj.UserData.MSTprojViewPlotObj.im_obj;
                        imObj.XData = [];
                        imObj.YData = [];
                        imObj.CData = [];
                    end
                else
                    switch upper(obj.objISV_proj.XY_COORDINATE_SYSTEM)
                        case {'PLANETOCENTRIC','LATLON'}
                            xy_coord = 'LATLON';
                        case 'NORTHEAST'
                            xy_coord = 'NE';
                        case 'IMAGEPIXELS'
                            xy_coord = 'PIXEL';
                        otherwise
                            
                    end
                    [selectMask_msldemc,srange,lrange] ...
                        = obj.MSTproj.mapper.get_msldemc_mask_base(x_demc,y_demc,xy_coord);
                    
                    if isempty(obj.objISVImage_ISV_proj_SelectMask)
                        obj.create_ax_ISV_proj_SelectMask();
                    end
                    
                    if isempty(cursor_obj.UserData.MSTprojViewPlotObj.im_obj)
                        col3 = reshape([1 0 0],[1,1,3]);
                        imObj = imagesc(obj.objISVImage_ISV_proj_SelectMask.ax,srange,lrange, ...
                            double(selectMask_msldemc).*col3,'AlphaData',double(selectMask_msldemc>0)*0.5);
                        obj.objISVImage_ISV_proj_SelectMask.imobj = [obj.objISVImage_ISV_proj_SelectMask.imobj imObj];
                        cursor_obj.UserData.MSTprojViewPlotObj.im_obj = ...
                            [cursor_obj.UserData.MSTprojViewPlotObj.im_obj imObj];
                    else
                        col3 = reshape([1 0 0],[1,1,3]);
                        imObj = cursor_obj.UserData.MSTprojViewPlotObj.im_obj;
                        imObj.XData = srange;
                        imObj.YData = lrange;
                        imObj.CData = double(selectMask_msldemc).*col3;
                        imObj.AlphaData = double(selectMask_msldemc>0)*0.5;
                     end
                    
                    
                end
            end
        end
        
        %%
        % =================================================================
        % ISV Projection view callback function
        % =================================================================
        function create_ax_MASTCAM_SelectMask(obj)
            obj.objISVImage_MASTCAM_SelectMask = obj.objENVIRasterMultview.obj_ISV.add_layer([],'name','Select Mask');
            obj.objISVImage_MASTCAM_SelectMask.ax.NextPlot = 'add';
        end
        function [out] = ISV_proj_BtnDwnFcn(obj,objAxes,eventData)
            [out] = obj.objISV_proj.image_BtnDwnFcn(objAxes,eventData);
            cursor_obj = out.cursor_obj;
            if isfield(cursor_obj.UserData,'MSTprojViewPlotObj')
            else
                cursor_obj.UserData.MSTprojViewPlotObj = MSTprojViewPlot();
            end
            obj.ISV_proj_update_cursor(cursor_obj);
            
        end
        
        function [out] = ISV_proj_WindowKeyPressFcn(obj,figobj,eventData)
            [out] = obj.objISV_proj.ISVWindowKeyPressFcn(figobj,eventData);
            cursor_obj = out.cursor_obj;
            switch eventData.Key
                case {'rightarrow','leftarrow','uparrow','downarrow'}
                    obj.ISV_proj_update_cursor(cursor_obj);
            end
        end
        
        function [] = ISV_proj_update_cursor(obj,cursor_obj)
            x_isv = cursor_obj.X; y_isv = cursor_obj.Y;
            switch upper(obj.objISV_proj.XY_COORDINATE_SYSTEM)
                case {'PLANETOCENTRIC','LATLON'}
                    x_demc = obj.MSTproj.mapper.msldemc_lon2x(x_isv);
                    y_demc = obj.MSTproj.mapper.msldemc_lat2y(y_isv);
                case 'NORTHEAST'
                    x_demc = obj.MSTproj.mapper.msldemc_easting2x(x_isv);
                    y_demc = obj.MSTproj.mapper.msldemc_northing2y(y_isv);
                case 'IMAGEPIXELS'
                    x_demc = x_isv;
                    y_demc = y_isv;
                otherwise

            end
            
            x_demc_int = round(x_demc); y_demc_int = round(y_demc);
            % if cursor is pointing the sky, then do not any further
            % processing
            [x_mst,y_mst] = obj.MSTproj.mapper.get_mastcam_index(x_demc_int,y_demc_int);
            if isempty(x_mst)
                if ~isempty(cursor_obj.UserData.MSTprojViewPlotObj.im_obj)
                    imObj = cursor_obj.UserData.MSTprojViewPlotObj.im_obj;
                    imObj.XData = [];
                    imObj.YData = [];
                    imObj.CData = [];
                end
            else
                [selectMask_mastcam,srange,lrange] = obj.MSTproj.mapper.get_mastcam_mask_base(x_mst,y_mst);

                if isempty(obj.objISVImage_MASTCAM_SelectMask)
                    obj.create_ax_MASTCAM_SelectMask();
                end

                if isempty(cursor_obj.UserData.MSTprojViewPlotObj.im_obj)
                    col3 = reshape([1 0 0],[1,1,3]);
                    imObj = imagesc(obj.objISVImage_MASTCAM_SelectMask.ax,srange,lrange, ...
                        double(selectMask_mastcam).*col3,'AlphaData',double(selectMask_mastcam>0)*0.5);
                    obj.objISVImage_MASTCAM_SelectMask.imobj = [obj.objISVImage_MASTCAM_SelectMask.imobj imObj];
                    cursor_obj.UserData.MSTprojViewPlotObj.im_obj = ...
                        [cursor_obj.UserData.MSTprojViewPlotObj.im_obj imObj];
                else
                    col3 = reshape([1 0 0],[1,1,3]);
                    imObj = cursor_obj.UserData.MSTprojViewPlotObj.im_obj;
                    imObj.XData = srange;
                    imObj.YData = lrange;
                    imObj.CData = double(selectMask_mastcam).*col3;
                    imObj.AlphaData = double(selectMask_mastcam>0)*0.5;
                end
                    

            end
        end
        
        %%
        
        
        
    end
end