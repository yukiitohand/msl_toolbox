classdef MASTCAMdata < HSI
    % MASTCAMdata class
    %  MASTCAM data
    
    properties
        lblpath
        lbl
        prop
        PRODUCT_ID
        RMC
        CAM_MDL
        ROVER_NAV
        CAM_MDL_GEO
        FILTER_NUMBER
        RADIANCE_FACTOR
        INSTRUMENT_ID
        L_im
        S_im
    end
    
    methods
        function obj = MASTCAMdata(basename,dirpath,varargin)
            
            rover_nav_ver = 'localized_interp';
            varargin_rmidx = [];
            if (rem(length(varargin),2)==1)
                error('Optional parameters should always go by pairs');
            else
                for i=1:2:(length(varargin)-1)
                    switch upper(varargin{i})
                        case {'ROVER_NAV_VERSION','ROVER_NAV_VER'}
                            rover_nav_ver = varargin{i+1};
                            varargin_rmidx = [varargin_rmidx i i+1];
                        otherwise
                            error('Unrecognized option: %s',varargin{i});
                    end
                end
            end
            varargin_HSI = varargin(setdiff(1:length(varargin),varargin_rmidx));
            obj@HSI(basename,dirpath,varargin_HSI{:});
            obj.lblpath = joinPath(dirpath,[basename '.lbl']);
            obj.lbl = pds3lblread(obj.lblpath);
            obj.hdr = mastcam_extract_imghdr_from_lbl(obj.lbl);
            obj.prop = getProp_basenameMASTCAM(basename);
            obj.CAM_MDL = get_cammera_model(obj);
            obj.RMC = get_rmc(obj);
            obj.PRODUCT_ID = obj.lbl.PRODUCT_ID;
            obj.ROVER_NAV  = obj.get_rover_nav('VERSION',rover_nav_ver);
            obj.get_filter_number();
            obj.get_instrument_id();
            obj.get_radiance_factor();
            obj.L_im = obj.hdr.lines;
            obj.S_im = obj.hdr.samples;
        end
        function [cahvor_mdl] = get_cammera_model(obj)
            [cahvor_mdl] = mastcam_get_cahvor_model(obj.lbl);
        end
        function [rmc_mdl] = get_rmc(obj)
            [rmc_mdl] = mastcam_get_RMC(obj.lbl);
        end
        function [rover_nav] = get_rover_nav(obj,varargin)
            [rover_nav,option] = msl_get_ROVER_NAV(obj,varargin{:});
            % [rover_nav] = get_ROVER_NAV_from_mslplc_view_csv(...
            %     obj.RMC.SITE,obj.RMC.DRIVE,obj.RMC.POSE,varargin{:});
        end
        
        function [] = get_CAM_MDL_GEO(obj)
            [imxy_direc_rov] = get_3d_pointing_from_CAHV_v2(...
                [obj.hdr.lines,obj.hdr.samples],obj.CAM_MDL);
            cmmdl_A_rov0 = obj.ROVER_NAV.rot_mat * obj.CAM_MDL.A';
            cmmdl_C_rov0 = obj.ROVER_NAV.rot_mat * obj.CAM_MDL.C';
            imxy_direc_rov_2d = reshape(imxy_direc_rov,[obj.hdr.lines*obj.hdr.samples,3])';
            imxy_direc_rov0_2d = obj.ROVER_NAV.rot_mat * imxy_direc_rov_2d;
            imxy_direc_rov0 = reshape(imxy_direc_rov0_2d',[obj.hdr.lines,obj.hdr.samples,3]);
            cmmdl_C_geo = cmmdl_C_rov0 + [obj.ROVER_NAV.NORTHING; 
                              obj.ROVER_NAV.EASTING;
                              -obj.ROVER_NAV.ELEVATION];
                          
            obj.CAM_MDL_GEO.A_rov0 = cmmdl_A_rov0;
            obj.CAM_MDL_GEO.C_rov0 = cmmdl_C_rov0;
            obj.CAM_MDL_GEO.C_geo  = cmmdl_C_geo;
            obj.CAM_MDL_GEO.imxy_direc_rov = imxy_direc_rov;
            obj.CAM_MDL_GEO.imxy_direc_rov0 = imxy_direc_rov0;
            
        end
        
        function get_instrument_id(obj)
            obj.INSTRUMENT_ID = obj.lbl.INSTRUMENT_ID;
        end
        
        function get_filter_number(obj)
            obj.FILTER_NUMBER = obj.lbl.GROUP_IMAGE_REQUEST_PARMS.FILTER_NUMBER;
        end
        
        function get_radiance_factor(obj)
            rf_raw = obj.lbl.GROUP_PROCESSING_PARMS.RADIANCE_SCALING_FACTOR;
            switch obj.INSTRUMENT_ID
                case 'MAST_LEFT'
                    switch obj.FILTER_NUMBER
                        case 0
                            obj.RADIANCE_FACTOR = rf_raw;
                        case 1
                            obj.RADIANCE_FACTOR = rf_raw(2);
                        case 2
                            obj.RADIANCE_FACTOR = rf_raw(3);
                        case 3
                            obj.RADIANCE_FACTOR = rf_raw(1);
                        case 4
                            obj.RADIANCE_FACTOR = rf_raw(1);
                        case 5
                            % should be same for all elements
                            obj.RADIANCE_FACTOR = rf_raw(1);
                        case 6
                            obj.RADIANCE_FACTOR = rf_raw(1);
                    end
                    
                case 'MAST_RIGHT'
                    switch obj.FILTER_NUMBER
                        case 0
                            obj.RADIANCE_FACTOR = rf_raw;
                        case 1
                            obj.RADIANCE_FACTOR = rf_raw(2);
                        case 2
                            obj.RADIANCE_FACTOR = rf_raw(3);
                        case 3
                            obj.RADIANCE_FACTOR = rf_raw(1);
                        case 4
                            % should be same for all elements
                            obj.RADIANCE_FACTOR = rf_raw(1);
                        case 5
                            % should be same for all elements
                            obj.RADIANCE_FACTOR = rf_raw(1);
                        case 6
                            % should be same for all elements
                            obj.RADIANCE_FACTOR = rf_raw(1);
                    end
            end
        end
        
        function [img_iof] = get_IoF(obj)
            if isempty(obj.img)
                obj.readimg();
            end
            switch obj.FILTER_NUMBER
                case 0
                    img_iof = reshape(obj.RADIANCE_FACTOR,[1,1,3]) .* obj.img;
                otherwise
                    img_iof = obj.RADIANCE_FACTOR .* obj.img;
            end
        end
        
        function delete(obj)
            delete(obj.CAM_MDL);
            delete(obj.CAM_MDL_GEO);
            delete(obj.ROVER_NAV);
            delete(obj.RMC);
        end   
    end
end

